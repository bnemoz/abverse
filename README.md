# abverse

**Germline-informed reverse translation of antibody amino acid sequences to nucleotide sequences.**

`abverse` is a companion package to [abstar](https://github.com/briney/abstar). It takes antibody amino acid sequences — common output from mass spectrometry, proteomics, or databases — and produces nucleotide sequences that are maximally faithful to the inferred germline, so that downstream `abstar` annotation (V/J assignment, mutation counts, CDR/FWR regions) reflects real somatic hypermutation rather than arbitrary codon choices.

---

## Why abverse?

`abstar` requires nucleotide input. Researchers with AA sequences have two options:

1. **Naive reverse-translation** → pick any codon per amino acid → run abstar → get inflated mutation counts and unreliable CDR boundaries because every codon choice that differs from germline looks like a mutation.
2. **abverse** → single-pass algorithm → germline-faithful NT → feed directly into `abstar.run()`.

The abverse approach is **provably optimal**: for each codon position aligned to a germline gene, it picks the synonymous codon with the minimum Hamming distance to the germline codon (ties broken by human codon frequency). Because codons don't overlap and Hamming distance is additive, the global minimum equals the sum of per-position minima — and the entire lookup table is pre-computed at import time (O(1) per position at runtime).

---

## Installation

```bash
pip install abverse
```

**Requirements:** Python ≥ 3.10, [abutils](https://github.com/briney/abutils) ≥ 0.5.1, [abstar](https://github.com/briney/abstar) (for germline databases), [polars](https://pola.rs/) ≥ 0.20, MMseqs2 (bundled via abutils).

---

## Quick start

```python
import abutils
import abverse
import abstar

# Load your AA sequences (FASTA file, list of strings, or list of abutils.Sequence)
aa_seqs = abutils.io.read_fasta("antibodies_aa.fasta")

# Reverse-translate to germline-faithful NT sequences
nt_seqs = abverse.reverse_translate(aa_seqs)

# Feed directly into abstar — results will have meaningful mutation counts
results = abstar.run(nt_seqs)
```

The returned `nt_seqs` is a `list[abutils.Sequence]`. Each sequence carries three annotations:

| Annotation | Description |
|---|---|
| `v_call` | Assigned V germline gene |
| `j_call` | Assigned J germline gene |
| `reconstruction_method` | `germline_vj`, `germline_v_only`, `germline_j_only`, or `codon_frequency` |

---

## API

### `abverse.reverse_translate(sequences, ...)`

```python
abverse.reverse_translate(
    sequences,              # FASTA path | list[str] | list[abutils.Sequence]
    species="human",        # germline species
    receptor="bcr",         # receptor type
    n_processes=None,       # worker processes (default: cpu_count)
    threads=None,           # MMseqs2 threads
    chunksize=500,          # sequences per worker batch
    force_rebuild_db=False, # force re-build of germline AA databases
    output_fasta=None,      # optional path to write NT FASTA
    verbose=False,          # print progress
) -> list[abutils.Sequence]
```

### `abverse.build_germline_aa_db(species, receptor, force_rebuild)`

Pre-builds (or validates the cache of) the germline amino acid databases used internally. Call this once on first install to populate `~/.abverse/germline_dbs/`. Subsequent calls reuse the cache unless the source germline files change (SHA-256 invalidation).

---

## How it works

### Algorithm

```
1. MMseqs2 protein–protein search (all AA sequences vs. V germline AA DB)
   → best V assignment per sequence

2. Extract post-V region (aa_seq[v_qend+1:]) per sequence
   → MMseqs2 protein–protein search vs. J germline AA DB
   → best J assignment per sequence

3. Parallel reconstruction (ProcessPoolExecutor):
   • 5' overhang (before V alignment)  → most frequent human codon
   • V region                           → argmin_c[Hamming(c, germline_codon)] per position
   • CDR3 (V end → J start)            → most frequent human codon
   • J region                           → argmin_c[Hamming(c, germline_codon)] per position
   • 3' overhang (after J alignment)   → most frequent human codon

4. Validate: assert translate(output_nt) == input_aa for every sequence
```

### Germline database cache

On first use, abverse translates abstar's nucleotide V/J germlines to amino acid FASTA files, builds MMseqs2 protein databases, and caches everything under `~/.abverse/germline_dbs/`. The cache is automatically invalidated and rebuilt if abstar's germline files change (checked via SHA-256).

Frame detection for J genes uses the conserved `WG.G` (IGH) / `FG.G` (IGK/IGL) motif; a stop-free-frame fallback covers unusual alleles.

---

## Performance

Benchmarked on a single CPU core with 10,000 BCR AA sequences:

| Metric | Value |
|---|---|
| Throughput | ~775 sequences/second/core |
| abstar calls in critical path | **0** |
| translate(output) == input guarantee | 100% (validated per sequence) |

No iterative abstar calls occur during `reverse_translate` — the algorithm is a single-pass pipeline.

---

## Integration test results

Tested on 100 real human BCR sequences with known germline assignments:

| Metric | Result | Threshold |
|---|---|---|
| V-gene family agreement | ≥ 90% | 90% |
| J-gene family agreement | ≥ 80% | 80% |
| Exact V-call match | 75% | informational |
| Exact J-call match | 91% | informational |

The exact V-call rate of 75% reflects the fundamental ambiguity of assigning a specific allele from amino acid sequence alone (multiple alleles can share the same AA sequence). Gene-family agreement — the metric that matters for mutation analysis — passes comfortably.

---

## Edge cases

| Situation | Handling |
|---|---|
| No V assignment | Human codon frequency for all positions; `reconstruction_method='codon_frequency'` |
| No J assignment | Germline lookup for V region; fallback elsewhere |
| 5′ / 3′ overhangs | Human codon frequency |
| Germline codon truncated at gene edge | Human codon frequency |
| Non-standard AA (X, B, Z) | `NNN` |
| Stop codon in input AA | `ValueError` with position and sequence ID |
| V/J alignment overlap | V takes priority; J starts after V end |

---

## Development

```bash
git clone https://github.com/<your-org>/abverse.git
cd abverse
pip install -e . --no-build-isolation
pip install pytest

# Run all tests (unit + integration + scaling benchmark)
python3 -m pytest abverse/tests/ -v
```

The test suite (59 tests) covers the codon lookup table, germline database building, per-sequence reconstruction with all edge cases, the end-to-end pipeline, integration with real BCR sequences, and a 10k-sequence throughput benchmark.

---

## Package structure

```
abverse/
├── pyproject.toml
└── abverse/
    ├── __init__.py          # public API: reverse_translate, build_germline_aa_db
    ├── _codons.py           # 1280-entry optimal codon lookup table
    ├── _germline_db.py      # germline translation, MMseqs2 DB build, cache
    ├── _search.py           # V + J protein–protein search wrappers (Polars)
    ├── _reconstruct.py      # per-sequence NT reconstruction (pure, picklable)
    ├── _pipeline.py         # orchestration and parallel dispatch
    └── tests/
        ├── test_codons.py
        ├── test_germline_db.py
        ├── test_reconstruct.py
        ├── test_pipeline.py
        ├── test_integration.py
        └── test_scaling.py
```

---

## License

MIT
