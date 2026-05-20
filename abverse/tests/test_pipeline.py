"""Smoke tests for _pipeline.py — orchestration layer."""

import os
import pytest
import abutils

import abverse


# Minimal real IGH sequence (IGHV3-23*01 / IGHJ4*02) — trimmed for speed
REAL_IgH_AA = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR"
    "WGQGTLVTVSS"
)


class TestReverseTranslate:
    def test_single_string_input(self):
        seqs = abverse.reverse_translate(["QVQLVQSGA"], n_processes=1, verbose=False)
        assert len(seqs) == 1
        assert isinstance(seqs[0], abutils.Sequence)

    def test_output_translates_back(self):
        aa = "EVQLVESGGGLVQ"
        seqs = abverse.reverse_translate([aa], n_processes=1)
        nt = str(seqs[0].sequence)
        assert len(nt) == len(aa) * 3
        # translate manually
        from abverse._reconstruct import _translate_nt
        assert _translate_nt(nt) == aa

    def test_fasta_input(self, tmp_path):
        fa = str(tmp_path / "test.fasta")
        with open(fa, "w") as fh:
            fh.write(">seq1\nQVQLVQSGA\n>seq2\nEVQLVESGGG\n")
        seqs = abverse.reverse_translate(fa, n_processes=1)
        assert len(seqs) == 2
        assert seqs[0].id == "seq1"
        assert seqs[1].id == "seq2"

    def test_real_igh_sequence(self):
        seqs = abverse.reverse_translate([REAL_IgH_AA], n_processes=1)
        assert len(seqs) == 1
        nt = str(seqs[0].sequence)
        assert len(nt) == len(REAL_IgH_AA) * 3
        from abverse._reconstruct import _translate_nt
        assert _translate_nt(nt) == REAL_IgH_AA

    def test_input_order_preserved(self):
        aas = ["QVQ", "MKP", "WGQ", "EVQ"]
        seqs = abverse.reverse_translate(aas, n_processes=1)
        from abverse._reconstruct import _translate_nt
        for aa, seq in zip(aas, seqs):
            assert _translate_nt(seq.sequence) == aa

    def test_output_fasta_written(self, tmp_path):
        out = str(tmp_path / "out.fasta")
        abverse.reverse_translate(["QVQLVQ"], output_fasta=out, n_processes=1)
        assert os.path.isfile(out)
        assert os.path.getsize(out) > 0

    def test_empty_input(self):
        seqs = abverse.reverse_translate([], n_processes=1)
        assert seqs == []
