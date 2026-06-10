# abverse — Loud failure on bad input

**Date:** 2026-06-10
**Status:** Approved (design)
**Author:** Benjamin Nemoz (with Claude Code)

## Problem

`reverse_translate` currently degrades silently when a sequence cannot be
reconstructed. Two distinct failure modes both end in the same trap:

1. **Invalid input residues.** Any amino acid outside the 20 standard residues
   (e.g. an `O` introduced by an OCR/transcription error) is not in the codon
   table. In `_reconstruct.py` only `X/B/Z/U` are special-cased to `NNN`;
   anything else falls through and produces `NNN` at that position too. The
   per-sequence validation `translate(output) == input` then fails because
   `translate("NNN") == "X"`, which is `!=` the original residue.

2. **Internal reconstruction errors.** When `reconstruct_sequence` raises (the
   validation `AssertionError`, or any unexpected exception), `_pipeline.py`
   catches it, substitutes `"N" * (3 * len(aa_seq))`, attaches a
   `reconstruction_error` annotation, and emits a `warnings.warn`.

The net effect: a single bad residue silently wipes the **entire** sequence to
`N`, and the only signal is a warning that is easily lost in a notebook. The
junk all-`N` output then propagates downstream and fails much later
(`abstar.run()` raised `MMseqs command failed ... no entry` on the empty
J-query FASTA), far from the real cause.

The promise in the README ("Non-standard AA (X, B, Z) → NNN, don't crash") is
also not actually upheld: because of the validation step, even the documented
`B/Z/U` codes trigger the whole-sequence wipe — only a literal `X` survives.

### Concrete trigger (motivating case)

`candidates.fasta` contained the residue `O` four times (HC pos 40; LC pos 4,
89, 105), each sitting where a `Q` belongs in a conserved antibody motif. Every
chain was silently returned as all-`N`; the failure only became visible as an
opaque MMseqs error inside `abstar.run()`.

## Goal

Make bad input **fail loudly and early**, with a single comprehensive error
that names every offending sequence — instead of silently producing degraded
output.

## Decisions

| Question | Decision |
|---|---|
| Failure timing | **Collect all** failures, then raise one comprehensive error — do not abort on the first bad sequence. |
| Allowed residues | **Only the 20 standard amino acids.** `X`, `B`, `Z`, `U`, `O`, `*`, and anything else are rejected. |
| Internal reconstruction errors | **Raise too.** Remove the silent all-`N` fallback entirely; collect and surface reconstruction exceptions. |

## Design

### 1. New public exception — `ReverseTranslationError`

Defined in a new module `abverse/_errors.py`, re-exported from
`abverse/__init__.py` (added to `__all__`).

```python
class ReverseTranslationError(Exception):
    """Raised when one or more sequences cannot be reverse-translated.

    Attributes
    ----------
    failures : list[dict]
        One entry per failed sequence, each with keys:
          - seq_id : str
          - kind   : "invalid_residue" | "reconstruction_error"
          - detail : str   (human-readable, e.g. "'O' at positions 4, 89, 105")
    """
    def __init__(self, failures: list[dict]):
        self.failures = failures
        super().__init__(self._render())

    def _render(self) -> str:
        ...
```

`_render()` produces a readable multi-line message, e.g.:

```
2 of 2 input sequences contain invalid residues
(only the 20 standard amino acids are accepted):
  MR-72_HC: 'O' at position 40
  MR-72_LC: 'O' at positions 4, 89, 105
```

For the reconstruction-error category the header reads
`N sequence(s) failed reconstruction:` followed by one line per `seq_id` and its
error detail.

The structured `.failures` list lets callers handle the error programmatically
rather than parsing the message.

### 2. Upfront input validation (the "collect all" gate)

A new pure helper in `abverse/_pipeline.py` (or a small validation helper
module):

```python
_STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

def _validate_residues(aa_seqs: dict[str, str]) -> list[dict]:
    """Return a list of invalid_residue failure dicts (empty if all valid)."""
```

For each sequence it records every position whose residue is not in
`_STANDARD_AA`, grouping by offending character for a compact message
(`'O' at positions 4, 89, 105`). `*` (stop) is reported here as just another
invalid residue.

It is called in `reverse_translate` **immediately after `_normalise_input`**,
so it covers every input path uniformly — FASTA, simple/AIRR/PairPlex CSV,
`abutils.Sequence`, and plain strings. If it returns any failures →
`raise ReverseTranslationError(failures)` **before** building the germline DB or
running any MMseqs2 search (no point spending work on invalid input). Because it
scans all sequences first, the single raised error names every bad sequence at
once, satisfying the "collect all" requirement.

### 3. Drop the silent all-`N` fallback

In both result-handling branches of `reverse_translate` (the single-process
path and the `ProcessPoolExecutor` path), reconstruction exceptions returned by
`_reconstruct_batch` are **collected** into a `reconstruction_error` failure
list instead of being converted to `"N" * (3 * len)` + `warnings.warn`. The
same applies to a whole-chunk future failure.

After all chunks are processed, if the reconstruction-failure list is non-empty
→ `raise ReverseTranslationError(failures)`.

Post-validation, these exceptions should essentially never fire for legitimate
input — so a raise here means a genuine abverse bug is surfacing immediately,
rather than emitting junk that fails downstream.

The `warnings` import and the `reconstruction_error` annotation path are
removed.

### 4. Simplify `_codons.py` / `_reconstruct.py`

Because input is guaranteed clean by the time reconstruction runs, the
non-standard-residue handling becomes dead code:

- `_reconstruct.py`: remove the `if aa in ("X", "B", "Z", "U"): nt_out[i] = "NNN"`
  branch in the per-position loop. The `"*" in aa_seq` guard at the top of
  `reconstruct_sequence` stays as belt-and-suspenders (cheap, defensive).
- `_codons.py`: `optimal_codon` and `fallback_codon` no longer special-case
  `X/B/Z/U` to `NNN`. `fallback_codon` raises a clear error (`KeyError` or
  `ValueError`) for an unknown residue instead of silently returning `NNN`, so
  any unexpected residue that reaches reconstruction is caught and surfaced as a
  `reconstruction_error` rather than producing silent junk.

`_translate_nt` (validation helper) is unchanged — it still maps unknown codons
to `"X"`, which is only used to compare against the input.

### 5. Documentation

- README: update the "Edge cases" table. The rows for `Non-standard AA (X, B, Z)`
  and `Stop codon in input AA` change to: rejected with a
  `ReverseTranslationError` naming the sequence(s) and offending positions.
  Note that only the 20 standard amino acids are accepted.
- `reverse_translate` docstring: replace the "Sequences that fail reconstruction
  are replaced with a Sequence containing 'N' …" paragraph with the new
  raise-on-failure contract, and document the `ReverseTranslationError` it may
  raise.

## Testing

- **Update** `test_reconstruct.py::test_nonstandard_aa_returns_nnn` → rename to
  `test_nonstandard_aa_raises` and assert reconstruction now raises (no `NNN`
  passthrough). `test_stop_codon_raises` stays (the defensive guard remains).
- **New, pipeline-level** (`test_pipeline.py`):
  - Input containing an invalid residue raises `ReverseTranslationError`.
  - The raised error's `.failures` / message names **all** bad sequences and
    their offending positions when multiple sequences are invalid (verifies
    "collect all").
  - Clean input is unaffected (existing happy-path tests stay green).
- **New** (`test_errors.py` or in `test_pipeline.py`): `ReverseTranslationError`
  renders a readable message and exposes a structured `.failures` list.
- The existing ~59-test suite otherwise remains green.

## Out of scope (YAGNI)

- No per-call `on_invalid` / `strict` flag — the strict, loud behavior is the
  single default. (Can be revisited later if a permissive mode is genuinely
  needed.)
- No automatic residue correction (e.g. `O → Q`). Fixing source data is the
  caller's responsibility.
- No changes to germline assignment, codon optimization, or the CSV format
  detection logic.

## Net effect

The motivating `O` case would raise, before any MMseqs2 work:

```
ReverseTranslationError: 2 of 2 input sequences contain invalid residues
(only the 20 standard amino acids are accepted):
  MR-72_HC: 'O' at position 40
  MR-72_LC: 'O' at positions 4, 89, 105
```

— pointing directly at the bad data instead of surfacing as an opaque MMseqs
failure inside `abstar.run()`.
