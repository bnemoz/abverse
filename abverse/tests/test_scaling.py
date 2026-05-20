"""Scaling benchmark: 10k sequences, measure throughput.

Generates a synthetic 10k dataset by repeating the 100-sequence test FASTA,
then verifies throughput and that no iterative abstar calls occur in the
critical path (reconstruction only, timing abverse.reverse_translate).

Run with:
    python3 -m pytest abverse/tests/test_scaling.py -v -s

Expected: ≥ 200 sequences/second on a single core (very conservative).
"""

import os
import time

import pytest
import abutils

import abverse
from abverse.tests.test_integration import _parse_paired_fasta

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
FASTA_PATH = os.path.join(DATA_DIR, "_100seqs_test.fasta")

pytestmark = pytest.mark.skipif(
    not os.path.isfile(FASTA_PATH),
    reason=f"Scaling benchmark FASTA not found: {FASTA_PATH}",
)

TARGET_N = 10_000
MIN_THROUGHPUT = 200  # sequences/second (conservative lower bound)


@pytest.fixture(scope="module")
def large_aa_seqs():
    """Repeat the 100-seq test set to reach TARGET_N sequences."""
    _, aa = _parse_paired_fasta(FASTA_PATH)
    base = list(aa.items())
    assert base, "No AA sequences parsed from test FASTA"

    seqs: list[abutils.Sequence] = []
    i = 0
    while len(seqs) < TARGET_N:
        sid, seq = base[i % len(base)]
        seqs.append(abutils.Sequence(seq, id=f"{sid}_{i // len(base)}"))
        i += 1
    return seqs[:TARGET_N]


class TestScaling:
    def test_throughput_single_core(self, large_aa_seqs):
        """reverse_translate 10k sequences on 1 process; assert ≥ MIN_THROUGHPUT seq/s."""
        n = len(large_aa_seqs)
        assert n == TARGET_N

        t0 = time.perf_counter()
        results = abverse.reverse_translate(large_aa_seqs, n_processes=1, verbose=False)
        elapsed = time.perf_counter() - t0

        throughput = n / elapsed
        print(f"\n10k throughput (1 process): {throughput:.0f} seq/s  ({elapsed:.1f}s total)")

        assert len(results) == n, f"Expected {n} results, got {len(results)}"
        assert throughput >= MIN_THROUGHPUT, (
            f"Throughput {throughput:.0f} seq/s < {MIN_THROUGHPUT} seq/s minimum"
        )

    def test_translation_roundtrip_10k(self, large_aa_seqs):
        """All 10k reconstructed sequences must translate back to their input AA."""
        from abverse._reconstruct import _translate_nt
        results = abverse.reverse_translate(large_aa_seqs, n_processes=1)
        aa_map = {s.id: str(s.sequence) for s in large_aa_seqs}
        errors: list[str] = []
        for seq in results:
            expected = aa_map[seq.id]
            got = _translate_nt(seq.sequence)
            if got != expected:
                errors.append(seq.id)
            if len(errors) > 5:
                break
        assert not errors, f"{len(errors)}+ sequences failed translation round-trip: {errors}"

    def test_no_abstar_in_critical_path(self, large_aa_seqs, monkeypatch):
        """Confirm abstar is never called inside reverse_translate."""
        import abstar as _abstar

        abstar_call_count = {"n": 0}
        original_run = _abstar.run

        def counting_run(*args, **kwargs):
            abstar_call_count["n"] += 1
            return original_run(*args, **kwargs)

        monkeypatch.setattr(_abstar, "run", counting_run)

        # Use a small subset for speed
        abverse.reverse_translate(large_aa_seqs[:50], n_processes=1)

        assert abstar_call_count["n"] == 0, (
            f"abstar.run() was called {abstar_call_count['n']} time(s) inside reverse_translate — "
            "iterative abstar calls violate the design constraint"
        )
