"""Smoke tests for _pipeline.py — orchestration layer."""

import os
import pytest
import abutils

import abverse
from abverse import ReverseTranslationError
from abverse._pipeline import _validate_residues


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


class TestValidateResidues:
    def test_clean_input_returns_no_failures(self):
        assert _validate_residues({"a": "QVQLVQSGA", "b": "EVQLVESGGG"}) == []

    def test_invalid_residue_reported_with_position(self):
        failures = _validate_residues({"MR-72_HC": "QVQLVQOGA"})
        assert len(failures) == 1
        f = failures[0]
        assert f["seq_id"] == "MR-72_HC"
        assert f["kind"] == "invalid_residue"
        assert "'O'" in f["detail"]
        # 'O' is the 7th residue (1-indexed)
        assert "position 7" in f["detail"]

    def test_multiple_positions_grouped_by_residue(self):
        # 'O' at 1-indexed positions 4, 9, 11
        failures = _validate_residues({"LC": "QVQOVQSGOGO"})
        assert len(failures) == 1
        assert failures[0]["detail"] == "'O' at positions 4, 9, 11"

    def test_stop_codon_is_invalid_residue(self):
        failures = _validate_residues({"x": "QVQ*VQ"})
        assert len(failures) == 1
        assert failures[0]["kind"] == "invalid_residue"
        assert "'*'" in failures[0]["detail"]

    def test_collects_all_bad_sequences(self):
        failures = _validate_residues({"hc": "QVQO", "lc": "EVBQ", "ok": "QVQL"})
        ids = {f["seq_id"] for f in failures}
        assert ids == {"hc", "lc"}


class TestReverseTranslateValidation:
    def test_invalid_residue_raises(self):
        with pytest.raises(ReverseTranslationError):
            abverse.reverse_translate(["QVQLVQOGA"], n_processes=1)

    def test_error_names_all_bad_sequences(self, tmp_path):
        fa = str(tmp_path / "bad.fasta")
        with open(fa, "w") as fh:
            fh.write(">MR-72_HC\nQVQLVQSGAO\n>MR-72_LC\nEVQOVESGGO\n")
        with pytest.raises(ReverseTranslationError) as exc_info:
            abverse.reverse_translate(fa, n_processes=1)
        err = exc_info.value
        ids = {f["seq_id"] for f in err.failures}
        assert ids == {"MR-72_HC", "MR-72_LC"}
        assert all(f["kind"] == "invalid_residue" for f in err.failures)

    def test_clean_input_unaffected(self):
        seqs = abverse.reverse_translate(["QVQLVQSGA"], n_processes=1)
        assert len(seqs) == 1

    def test_reconstruction_error_raises_not_silent_n(self, monkeypatch):
        # Reconstruction exceptions must be collected and surfaced as a
        # ReverseTranslationError — never silently replaced with all-N output.
        import abverse._pipeline as pipeline

        def fake_batch(records, *args, **kwargs):
            return [ValueError("boom") for _ in records]

        monkeypatch.setattr(pipeline, "_reconstruct_batch", fake_batch)

        with pytest.raises(ReverseTranslationError) as exc_info:
            abverse.reverse_translate(["QVQLVQSGA"], n_processes=1)
        failures = exc_info.value.failures
        assert len(failures) == 1
        assert failures[0]["kind"] == "reconstruction_error"
        assert "boom" in failures[0]["detail"]
