"""Integration test: reverse_translate → abstar.run() → verify V/J calls match.

Uses the 100-sequence test FASTA in tests/data/_100seqs_test.fasta.
Each antibody appears twice in the file (same ID): first as NT, then as AA.
Ground-truth V/J calls are obtained by running abstar on the original NT sequences.
Pass criterion: ≥ 90% V-gene family agreement and ≥ 80% J-gene family agreement.
"""

import os
import re
from collections import defaultdict
from typing import Optional

import pytest
import abutils
import abstar

import abverse

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
FASTA_PATH = os.path.join(DATA_DIR, "_100seqs_test.fasta")

# Skip entire module if the test FASTA is absent
pytestmark = pytest.mark.skipif(
    not os.path.isfile(FASTA_PATH),
    reason=f"Integration test FASTA not found: {FASTA_PATH}",
)


# ── FASTA parser ──────────────────────────────────────────────────────────────

_NT_CHARS = set("ACGTacgt")
_AA_CHARS = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")


def _is_nt(seq: str) -> bool:
    return all(c in _NT_CHARS for c in seq)


def _parse_paired_fasta(path: str) -> tuple[dict[str, str], dict[str, str]]:
    """Parse a FASTA where each sequence appears twice (NT then AA) with the same ID.

    Returns (nt_seqs, aa_seqs) as {id: sequence} dicts.
    Duplicate IDs are resolved by inspecting sequence content.
    """
    nt: dict[str, str] = {}
    aa: dict[str, str] = {}
    seen: dict[str, int] = defaultdict(int)

    with open(path) as fh:
        current_id: Optional[str] = None
        current_parts: list[str] = []

        def _store(sid: str, parts: list[str]) -> None:
            seq = "".join(parts).strip()
            if not seq:
                return
            seen[sid] += 1
            if _is_nt(seq):
                nt[sid] = seq
            else:
                # AA sequence — deduplicate key if NT was already stored under same id
                key = sid if sid not in aa else f"{sid}_{seen[sid]}"
                aa[sid] = seq

        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    _store(current_id, current_parts)
                current_id = line[1:].split()[0]
                current_parts = []
            elif line:
                current_parts.append(line)
        if current_id is not None:
            _store(current_id, current_parts)

    return nt, aa


# ── Gene family extraction ────────────────────────────────────────────────────

def _gene_family(call: Optional[str]) -> Optional[str]:
    """Return the V/J gene family, e.g. 'IGHV3' from 'IGHV3-23*01__homo_sapiens'."""
    if not call:
        return None
    # Strip species suffix if present
    call = call.split("__")[0]
    # Match IGHV3, IGKV1, IGLJ2, etc.
    m = re.match(r"(IG[HKL][VDJ]\d+)", call, re.IGNORECASE)
    return m.group(1).upper() if m else call.split("-")[0].upper()


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def parsed_seqs():
    return _parse_paired_fasta(FASTA_PATH)


@pytest.fixture(scope="module")
def nt_seqs(parsed_seqs):
    nt, _ = parsed_seqs
    return nt


@pytest.fixture(scope="module")
def aa_seqs(parsed_seqs):
    _, aa = parsed_seqs
    return aa


@pytest.fixture(scope="module")
def ground_truth_calls(nt_seqs):
    """Run abstar on original NT sequences; return {seq_id: (v_call, j_call)}."""
    seqs = [abutils.Sequence(seq, id=sid) for sid, seq in nt_seqs.items()]
    results = abstar.run(seqs)
    if isinstance(results, abutils.Sequence):
        results = [results]
    calls: dict[str, tuple[Optional[str], Optional[str]]] = {}
    for result in results:
        sid = result["sequence_id"]
        calls[sid] = (result.get("v_call"), result.get("j_call"))
    return calls


@pytest.fixture(scope="module")
def revtrans_results(aa_seqs):
    """Run abverse.reverse_translate on the 100 AA sequences."""
    aa_list = [abutils.Sequence(seq, id=sid) for sid, seq in aa_seqs.items()]
    return abverse.reverse_translate(aa_list, n_processes=1)


@pytest.fixture(scope="module")
def revtrans_calls(revtrans_results):
    """Run abstar on the reverse-translated NT sequences; return {seq_id: (v_call, j_call)}."""
    results = abstar.run(revtrans_results)
    if isinstance(results, abutils.Sequence):
        results = [results]
    calls: dict[str, tuple[Optional[str], Optional[str]]] = {}
    for result in results:
        sid = result["sequence_id"]
        calls[sid] = (result.get("v_call"), result.get("j_call"))
    return calls


# ── Tests ─────────────────────────────────────────────────────────────────────

class TestIntegration:
    def test_fasta_parsed(self, nt_seqs, aa_seqs):
        assert len(nt_seqs) >= 50, f"Expected ≥50 NT seqs, got {len(nt_seqs)}"
        assert len(aa_seqs) >= 50, f"Expected ≥50 AA seqs, got {len(aa_seqs)}"

    def test_all_aa_seqs_produce_nt(self, revtrans_results, aa_seqs):
        assert len(revtrans_results) == len(aa_seqs)
        for seq in revtrans_results:
            assert len(seq.sequence) % 3 == 0
            assert len(seq.sequence) > 0

    def test_reverse_translate_output_translates_back(self, revtrans_results, aa_seqs):
        from abverse._reconstruct import _translate_nt
        errors: list[str] = []
        for seq in revtrans_results:
            expected_aa = aa_seqs[seq.id]
            got = _translate_nt(seq.sequence)
            if got != expected_aa:
                errors.append(
                    f"{seq.id}: translate(output) != input\n"
                    f"  expected : {expected_aa}\n"
                    f"  got      : {got}"
                )
        assert not errors, f"{len(errors)} sequences failed translation round-trip:\n" + "\n".join(errors[:3])

    def test_v_gene_family_agreement(self, ground_truth_calls, revtrans_calls):
        """≥ 90% of sequences should have matching V-gene family after round-trip."""
        shared = set(ground_truth_calls) & set(revtrans_calls)
        assert shared, "No shared sequence IDs between ground-truth and revtrans calls"

        matches = 0
        total = 0
        mismatches: list[str] = []
        for sid in shared:
            gt_v, _ = ground_truth_calls[sid]
            rt_v, _ = revtrans_calls[sid]
            if gt_v is None:
                continue  # skip sequences abstar couldn't annotate
            total += 1
            gt_fam = _gene_family(gt_v)
            rt_fam = _gene_family(rt_v)
            if gt_fam == rt_fam:
                matches += 1
            else:
                mismatches.append(f"  {sid}: GT={gt_v!r}  RT={rt_v!r}")

        rate = matches / total if total else 0
        assert rate >= 0.90, (
            f"V-gene family agreement {rate:.1%} < 90% ({matches}/{total} matched)\n"
            + "\n".join(mismatches[:10])
        )

    def test_j_gene_family_agreement(self, ground_truth_calls, revtrans_calls):
        """≥ 80% of sequences should have matching J-gene family after round-trip."""
        shared = set(ground_truth_calls) & set(revtrans_calls)
        matches = 0
        total = 0
        mismatches: list[str] = []
        for sid in shared:
            _, gt_j = ground_truth_calls[sid]
            _, rt_j = revtrans_calls[sid]
            if gt_j is None:
                continue
            total += 1
            gt_fam = _gene_family(gt_j)
            rt_fam = _gene_family(rt_j)
            if gt_fam == rt_fam:
                matches += 1
            else:
                mismatches.append(f"  {sid}: GT={gt_j!r}  RT={rt_j!r}")

        rate = matches / total if total else 0
        assert rate >= 0.80, (
            f"J-gene family agreement {rate:.1%} < 80% ({matches}/{total} matched)\n"
            + "\n".join(mismatches[:10])
        )

    def test_exact_v_call_agreement(self, ground_truth_calls, revtrans_calls):
        """Report (don't assert) exact V-call agreement rate — informational."""
        shared = set(ground_truth_calls) & set(revtrans_calls)
        matches = sum(
            1 for sid in shared
            if ground_truth_calls[sid][0] and
               ground_truth_calls[sid][0] == revtrans_calls[sid][0]
        )
        total = sum(1 for sid in shared if ground_truth_calls[sid][0])
        rate = matches / total if total else 0
        print(f"\nExact V-call agreement: {matches}/{total} = {rate:.1%}")
        # No assertion — this is an informational metric

    def test_exact_j_call_agreement(self, ground_truth_calls, revtrans_calls):
        """Report (don't assert) exact J-call agreement rate — informational."""
        shared = set(ground_truth_calls) & set(revtrans_calls)
        matches = sum(
            1 for sid in shared
            if ground_truth_calls[sid][1] and
               ground_truth_calls[sid][1] == revtrans_calls[sid][1]
        )
        total = sum(1 for sid in shared if ground_truth_calls[sid][1])
        rate = matches / total if total else 0
        print(f"\nExact J-call agreement: {matches}/{total} = {rate:.1%}")
