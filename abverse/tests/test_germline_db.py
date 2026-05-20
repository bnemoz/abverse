"""Tests for _germline_db.py — germline translation and database building."""

import os

import pytest

from abverse._germline_db import (
    _build_j_aa_fasta,
    _build_v_aa_fasta,
    _detect_j_frame,
    _read_nt_fasta,
    _translate_nt,
    build_germline_aa_db,
)


@pytest.fixture(scope="module")
def v_fasta_path():
    import abstar
    base = os.path.dirname(abstar.__file__)
    return os.path.join(base, "germline_dbs", "bcr", "human", "ungapped", "v.fasta")


@pytest.fixture(scope="module")
def j_fasta_path():
    import abstar
    base = os.path.dirname(abstar.__file__)
    return os.path.join(base, "germline_dbs", "bcr", "human", "ungapped", "j.fasta")


@pytest.fixture(scope="module")
def v_nt_seqs(v_fasta_path):
    return _read_nt_fasta(v_fasta_path)


@pytest.fixture(scope="module")
def j_nt_seqs(j_fasta_path):
    return _read_nt_fasta(j_fasta_path)


class TestReadFasta:
    def test_v_genes_loaded(self, v_nt_seqs):
        assert len(v_nt_seqs) > 0

    def test_j_genes_loaded(self, j_nt_seqs):
        assert len(j_nt_seqs) > 0

    def test_sequences_are_nt(self, v_nt_seqs):
        for name, seq in list(v_nt_seqs.items())[:5]:
            assert all(c in "ACGTacgt" for c in seq), f"Non-ACGT in {name}"


class TestTranslateNt:
    def test_methionine(self):
        assert _translate_nt("ATG") == "M"

    def test_stop_codon(self):
        assert _translate_nt("TAA") == "*"

    def test_full_codon(self):
        assert _translate_nt("ATGAAACCC") == "MKP"


class TestDetectJFrame:
    def test_igh_j_gene(self, j_nt_seqs):
        # IGH J4*02: ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG
        # Frame 1: TTLTLAGGTP… no WG/FG motif
        # Need to find a gene with WG.G motif
        for name, nt in j_nt_seqs.items():
            frame = _detect_j_frame(nt)
            assert frame in (1, 2, 3), f"Invalid frame {frame} for {name}"

    def test_all_j_genes_have_stop_free_frame(self, j_nt_seqs):
        for name, nt in j_nt_seqs.items():
            frame = _detect_j_frame(nt)
            offset = frame - 1
            aa = _translate_nt(nt[offset:])
            # The translation up to the first stop should be at least 5 AA
            clean = aa.split("*")[0]
            assert len(clean) >= 5, f"J gene {name} (frame {frame}) too short: {aa!r}"


class TestBuildVaaFasta:
    def test_all_v_genes_translate_clean(self, v_nt_seqs, tmp_path):
        out = str(tmp_path / "v_aa.fasta")
        frame_map = _build_v_aa_fasta(v_nt_seqs, out)
        # Read back and verify no stop codons in the main body
        seqs = _read_nt_fasta(out)
        for name, aa in seqs.items():
            assert "*" not in aa, f"Stop codon in V gene AA: {name}: {aa!r}"
        assert len(frame_map) == len(v_nt_seqs)
        assert all(f == 1 for f in frame_map.values())

    def test_fasta_file_written(self, v_nt_seqs, tmp_path):
        out = str(tmp_path / "v_aa.fasta")
        _build_v_aa_fasta(v_nt_seqs, out)
        assert os.path.isfile(out)
        assert os.path.getsize(out) > 0


class TestBuildJaaFasta:
    def test_all_j_genes_translate_clean(self, j_nt_seqs, tmp_path):
        out = str(tmp_path / "j_aa.fasta")
        frame_map = _build_j_aa_fasta(j_nt_seqs, out)
        seqs = _read_nt_fasta(out)
        for name, aa in seqs.items():
            assert "*" not in aa, f"Stop codon in J gene AA: {name}: {aa!r}"
        assert len(frame_map) == len(j_nt_seqs)


class TestBuildGermlineAaDb:
    def test_builds_and_caches(self, tmp_path, monkeypatch):
        # Patch cache root to tmp_path
        import abverse._germline_db as gdb
        monkeypatch.setattr(gdb, "_CACHE_ROOT", str(tmp_path))
        db = build_germline_aa_db(species="human", receptor="bcr")
        assert "v_db_path" in db
        assert "j_db_path" in db
        assert "v_nt_seqs" in db
        assert "j_nt_seqs" in db
        assert "v_frame_map" in db
        assert "j_frame_map" in db
        assert os.path.isfile(db["v_db_path"])
        assert os.path.isfile(db["j_db_path"])

    def test_cache_reuse(self, tmp_path, monkeypatch):
        import abverse._germline_db as gdb
        monkeypatch.setattr(gdb, "_CACHE_ROOT", str(tmp_path))
        db1 = build_germline_aa_db(species="human", receptor="bcr")
        db2 = build_germline_aa_db(species="human", receptor="bcr")
        assert db1["v_db_path"] == db2["v_db_path"]

    def test_force_rebuild(self, tmp_path, monkeypatch):
        import abverse._germline_db as gdb
        monkeypatch.setattr(gdb, "_CACHE_ROOT", str(tmp_path))
        db1 = build_germline_aa_db(species="human", receptor="bcr")
        db2 = build_germline_aa_db(species="human", receptor="bcr", force_rebuild=True)
        assert db2["v_db_path"] == db1["v_db_path"]
