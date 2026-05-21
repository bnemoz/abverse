"""Tests for _csv.py — CSV parsing for simple, AIRR, and PairPlex formats."""

import os
import textwrap

import pytest

from abverse._csv import (
    _extract_airr_coords,
    _parse_simple,
    detect_csv_format,
    parse_csv,
)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _write(tmp_path, name: str, content: str) -> str:
    p = tmp_path / name
    p.write_text(textwrap.dedent(content))
    return str(p)


# ── detect_csv_format ─────────────────────────────────────────────────────────

class TestDetectCsvFormat:
    def test_simple(self):
        assert detect_csv_format(["id", "sequence"]) == "simple"

    def test_airr_with_sequence_aa(self):
        cols = ["sequence_id", "sequence_aa", "v_call", "j_call"]
        assert detect_csv_format(cols) == "airr"

    def test_airr_with_sequence_only(self):
        cols = ["sequence_id", "sequence", "v_call"]
        assert detect_csv_format(cols) == "airr"

    def test_pairplex_suffix_0(self):
        cols = ["name", "sequence_aa:0", "sequence_aa:1", "v_call:0"]
        assert detect_csv_format(cols) == "pairplex"

    def test_pairplex_suffix_1_only(self):
        # Even a single :1 column triggers pairplex detection
        cols = ["name", "sequence_aa:1"]
        assert detect_csv_format(cols) == "pairplex"

    def test_ambiguous_no_sequence_id(self):
        # Has v_call but no sequence_id → not AIRR → simple
        assert detect_csv_format(["v_call", "sequence_aa"]) == "simple"


# ── _extract_airr_coords ──────────────────────────────────────────────────────

class TestExtractAirrCoords:
    def test_full_vj(self):
        row = {
            "v_sequence_start": 1,  "v_sequence_end": 285,
            "v_germline_start": 1,
            "j_sequence_start": 310, "j_sequence_end": 345,
            "j_germline_start": 1,
        }
        c = _extract_airr_coords(row)
        assert c["v_qstart"] == 0
        assert c["v_qend"]   == 94   # (285-1)//3
        assert c["v_tstart"] == 0
        # j relative to post-V: abs_start = (310-1)//3 = 103
        # v_qend = 94 → post-V starts at 95 → j_qstart = 103 - 95 = 8
        assert c["j_qstart"] == 8
        assert c["j_qend"]   == (345 - 1) // 3 - 95
        assert c["j_tstart"] == 0

    def test_v_only(self):
        row = {
            "v_sequence_start": 1, "v_sequence_end": 285,
            "v_germline_start": 1,
            "j_sequence_start": None, "j_sequence_end": None, "j_germline_start": None,
        }
        c = _extract_airr_coords(row)
        assert "v_qstart" in c
        assert "j_qstart" not in c

    def test_j_only(self):
        row = {
            "v_sequence_start": None, "v_sequence_end": None, "v_germline_start": None,
            "j_sequence_start": 310, "j_sequence_end": 345,
            "j_germline_start": 4,
        }
        c = _extract_airr_coords(row)
        assert "v_qstart" not in c
        # No V → j_qstart is absolute (j_abs_start when no v_qend)
        assert c["j_qstart"] == (310 - 1) // 3
        assert c["j_tstart"] == (4 - 1) // 3

    def test_all_none(self):
        row = {k: None for k in [
            "v_sequence_start", "v_sequence_end", "v_germline_start",
            "j_sequence_start", "j_sequence_end", "j_germline_start",
        ]}
        assert _extract_airr_coords(row) == {}

    def test_non_unit_v_germline_start(self):
        row = {
            "v_sequence_start": 1,  "v_sequence_end": 60,
            "v_germline_start": 10,   # germline starts mid-gene
            "j_sequence_start": None, "j_sequence_end": None, "j_germline_start": None,
        }
        c = _extract_airr_coords(row)
        assert c["v_tstart"] == (10 - 1) // 3  # == 3


# ── Simple CSV ────────────────────────────────────────────────────────────────

class TestSimpleCsv:
    def test_basic(self, tmp_path):
        p = _write(tmp_path, "simple.csv", """\
            id,sequence
            seq1,EVQLVESGG
            seq2,QVQLVQSGA
        """)
        aa, ann = parse_csv(p)
        assert aa == {"seq1": "EVQLVESGG", "seq2": "QVQLVQSGA"}
        assert ann == {}

    def test_no_header_still_works(self, tmp_path):
        # When header names are generic (non-AIRR), treated as simple CSV
        p = _write(tmp_path, "noh.csv", """\
            name,seq
            abc,MKPWST
        """)
        aa, _ = parse_csv(p)
        assert "abc" in aa
        assert aa["abc"] == "MKPWST"

    def test_lowercase_sequence_uppercased(self, tmp_path):
        p = _write(tmp_path, "lower.csv", """\
            id,seq
            s1,evqlves
        """)
        aa, _ = parse_csv(p)
        assert aa["s1"] == "EVQLVES"

    def test_empty_rows_skipped(self, tmp_path):
        p = _write(tmp_path, "empty.csv", """\
            id,seq
            s1,EVQLVES
            ,
            s2,QVQLVQ
        """)
        aa, _ = parse_csv(p)
        assert len(aa) == 2
        assert "s1" in aa and "s2" in aa

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_csv(str(tmp_path / "nonexistent.csv"))

    def test_tsv_extension_parsed_correctly(self, tmp_path):
        p = _write(tmp_path, "data.tsv", """\
            id\tsequence
            s1\tMKPWST
        """)
        aa, _ = parse_csv(p)
        assert aa.get("s1") == "MKPWST"


# ── AIRR CSV ──────────────────────────────────────────────────────────────────

class TestAirrCsv:
    def _make_airr_csv(self, tmp_path, rows: str, extra_header: str = "") -> str:
        header = f"sequence_id,sequence_aa,v_call,j_call{extra_header}"
        return _write(tmp_path, "airr.csv", f"{header}\n{rows}")

    def test_basic_gene_calls(self, tmp_path):
        p = self._make_airr_csv(
            tmp_path,
            "s1,EVQLVESGG,IGHV3-23*01,IGHJ4*02\n"
            "s2,QVQLVQSGA,IGHV1-2*02,IGHJ6*02",
        )
        aa, ann = parse_csv(p)
        assert aa == {"s1": "EVQLVESGG", "s2": "QVQLVQSGA"}
        assert ann["s1"]["v_call"] == "IGHV3-23*01"
        assert ann["s1"]["j_call"] == "IGHJ4*02"
        assert ann["s2"]["v_call"] == "IGHV1-2*02"

    def test_gene_calls_with_coordinates(self, tmp_path):
        header_extra = (
            ",v_sequence_start,v_sequence_end,v_germline_start"
            ",j_sequence_start,j_sequence_end,j_germline_start"
        )
        row = "s1,EVQLVESGG,IGHV3-23*01,IGHJ4*02,1,27,1,31,45,1"
        p = self._make_airr_csv(tmp_path, row, extra_header=header_extra)
        _, ann = parse_csv(p)
        assert "v_qstart" in ann["s1"]
        assert ann["s1"]["v_qstart"] == 0
        assert ann["s1"]["v_qend"]   == (27 - 1) // 3
        assert ann["s1"]["v_tstart"] == 0
        assert "j_qstart" in ann["s1"]

    def test_no_gene_call_cols_gives_no_annotations(self, tmp_path):
        p = _write(tmp_path, "airr_noc.csv", """\
            sequence_id,sequence_aa
            s1,EVQLVESGG
        """)
        _, ann = parse_csv(p)
        assert ann == {}

    def test_missing_sequence_aa_raises(self, tmp_path):
        p = _write(tmp_path, "airr_bad.csv", """\
            sequence_id,sequence,v_call
            s1,ATGCCC,IGHV1-2*02
        """)
        with pytest.raises(ValueError, match="sequence_aa"):
            parse_csv(p)

    def test_null_gene_calls_stored_as_none(self, tmp_path):
        p = self._make_airr_csv(
            tmp_path,
            "s1,EVQLVESGG,,\n"
        )
        _, ann = parse_csv(p)
        assert ann["s1"]["v_call"] is None
        assert ann["s1"]["j_call"] is None

    def test_imgt_gaps_stripped(self, tmp_path):
        p = self._make_airr_csv(tmp_path, "s1,EVQ.LVE..SGG,IGHV3-23*01,IGHJ4*02")
        aa, _ = parse_csv(p)
        assert aa["s1"] == "EVQLVESGG"

    def test_tsv_airr(self, tmp_path):
        p = _write(tmp_path, "airr.tsv", """\
            sequence_id\tsequence_aa\tv_call\tj_call
            s1\tEVQLVESGG\tIGHV3-23*01\tIGHJ4*02
        """)
        aa, ann = parse_csv(p)
        assert aa["s1"] == "EVQLVESGG"
        assert ann["s1"]["v_call"] == "IGHV3-23*01"

    def test_multiple_sequences_all_annotated(self, tmp_path):
        rows = "\n".join(
            f"s{i},{'EVQLVESGG'[:i+3]},IGHV3-23*01,IGHJ4*02" for i in range(5)
        )
        p = self._make_airr_csv(tmp_path, rows)
        aa, ann = parse_csv(p)
        assert len(aa) == 5
        assert all(d["v_call"] == "IGHV3-23*01" for d in ann.values())


# ── PairPlex CSV ──────────────────────────────────────────────────────────────

class TestPairplexCsv:
    def _make_pairplex_csv(self, tmp_path, rows: str, extra_header: str = "") -> str:
        header = (
            f"name,sequence_id:0,sequence_aa:0,v_call:0,j_call:0,"
            f"locus:0,sequence_id:1,sequence_aa:1,v_call:1,j_call:1,locus:1"
            f"{extra_header}"
        )
        return _write(tmp_path, "pairplex.csv", f"{header}\n{rows}")

    def test_basic_two_chains(self, tmp_path):
        row = (
            "pair1,pair1_heavy,EVQLVESGG,IGHV3-23*01,IGHJ4*02,IGH,"
            "pair1_light,QSVLTQPPA,IGLV2-14*01,IGLJ2*01,IGL"
        )
        p = self._make_pairplex_csv(tmp_path, row)
        aa, ann = parse_csv(p)
        assert "pair1_heavy" in aa
        assert "pair1_light" in aa
        assert aa["pair1_heavy"] == "EVQLVESGG"
        assert aa["pair1_light"] == "QSVLTQPPA"
        assert ann["pair1_heavy"]["v_call"] == "IGHV3-23*01"
        assert ann["pair1_light"]["v_call"] == "IGLV2-14*01"

    def test_chain_label_from_locus_when_no_sequence_id(self, tmp_path):
        # No sequence_id:0 / :1 columns → fall back to name + chain label
        header = "name,sequence_aa:0,v_call:0,j_call:0,locus:0,sequence_aa:1,v_call:1,j_call:1,locus:1"
        row = "ab1,EVQLVESGG,IGHV3-23*01,IGHJ4*02,IGH,QSVLTQPPA,IGLV2-14*01,IGLJ2*01,IGL"
        p = _write(tmp_path, "pairplex2.csv", f"{header}\n{row}")
        aa, _ = parse_csv(p)
        assert "ab1_H" in aa
        assert "ab1_L" in aa

    def test_chain_label_fallback_no_locus(self, tmp_path):
        header = "name,sequence_aa:0,v_call:0,j_call:0,sequence_aa:1,v_call:1,j_call:1"
        row = "ab1,EVQLVESGG,IGHV3-23*01,IGHJ4*02,QSVLTQPPA,IGLV2-14*01,IGLJ2*01"
        p = _write(tmp_path, "pairplex3.csv", f"{header}\n{row}")
        aa, _ = parse_csv(p)
        assert "ab1_chain0" in aa
        assert "ab1_chain1" in aa

    def test_pairplex_with_coordinates(self, tmp_path):
        coord_header = (
            ",v_sequence_start:0,v_sequence_end:0,v_germline_start:0"
            ",j_sequence_start:0,j_sequence_end:0,j_germline_start:0"
            ",v_sequence_start:1,v_sequence_end:1,v_germline_start:1"
            ",j_sequence_start:1,j_sequence_end:1,j_germline_start:1"
        )
        row = (
            "pair1,h1,EVQLVESGG,IGHV3-23*01,IGHJ4*02,IGH,"
            "l1,QSVLTQPPA,IGLV2-14*01,IGLJ2*01,IGL,"
            "1,27,1,31,45,1,"
            "1,27,1,31,45,1"
        )
        p = self._make_pairplex_csv(tmp_path, row, extra_header=coord_header)
        _, ann = parse_csv(p)
        assert "v_qstart" in ann["h1"]
        assert "v_qstart" in ann["l1"]
        assert ann["h1"]["v_qstart"] == 0

    def test_multiple_pairs(self, tmp_path):
        rows = "\n".join(
            f"pair{i},h{i},EVQLVES,IGHV3-23*01,IGHJ4*02,IGH,"
            f"l{i},QSVLTQ,IGLV2-14*01,IGLJ2*01,IGL"
            for i in range(3)
        )
        p = self._make_pairplex_csv(tmp_path, rows)
        aa, ann = parse_csv(p)
        assert len(aa) == 6  # 3 pairs × 2 chains
        assert all(d.get("v_call") is not None for d in ann.values())

    def test_missing_chain_1_skipped(self, tmp_path):
        # Only :0 columns present
        header = "name,sequence_aa:0,v_call:0,j_call:0,locus:0"
        row = "pair1,EVQLVESGG,IGHV3-23*01,IGHJ4*02,IGH"
        p = _write(tmp_path, "pairplex_half.csv", f"{header}\n{row}")
        aa, _ = parse_csv(p)
        assert len(aa) == 1  # only heavy chain

    def test_null_sequence_aa_skipped(self, tmp_path):
        row = (
            "pair1,h1,,IGHV3-23*01,IGHJ4*02,IGH,"
            "l1,QSVLTQPPA,IGLV2-14*01,IGLJ2*01,IGL"
        )
        p = self._make_pairplex_csv(tmp_path, row)
        aa, _ = parse_csv(p)
        # heavy chain has null sequence_aa → skipped
        assert "h1" not in aa
        assert "l1" in aa


# ── Pipeline integration ──────────────────────────────────────────────────────

class TestPipelineCsvInput:
    """Verify that reverse_translate() accepts CSV paths end-to-end."""

    def test_simple_csv_round_trip(self, tmp_path):
        import abverse
        from abverse._reconstruct import _translate_nt

        p = _write(tmp_path, "simple.csv", """\
            id,sequence
            seq1,EVQLVESGG
            seq2,QVQLVQ
        """)
        results = abverse.reverse_translate(p, n_processes=1)
        assert len(results) == 2
        ids = {r.id for r in results}
        assert "seq1" in ids
        assert "seq2" in ids
        for r in results:
            aa_in = "EVQLVESGG" if r.id == "seq1" else "QVQLVQ"
            assert _translate_nt(r.sequence) == aa_in

    def test_airr_csv_skips_mmseqs2_for_annotated(self, tmp_path):
        """With pre-annotated AIRR input, MMseqs2 is bypassed: no MMseqs2 binary needed."""
        import abverse
        from abverse._reconstruct import _translate_nt
        from unittest.mock import patch

        p = _write(tmp_path, "airr.csv", """\
            sequence_id,sequence_aa,v_call,j_call
            s1,EVQLVESGG,IGHV3-23*01,IGHJ4*02
        """)

        # Patch search functions to assert they are NOT called
        with patch("abverse._pipeline.search_v_germline") as mock_v, \
             patch("abverse._pipeline.search_j_germline") as mock_j:
            results = abverse.reverse_translate(p, n_processes=1)
            mock_v.assert_not_called()
            mock_j.assert_not_called()

        assert len(results) == 1
        assert _translate_nt(results[0].sequence) == "EVQLVESGG"

    def test_airr_csv_gene_calls_preserved(self, tmp_path):
        import abverse

        p = _write(tmp_path, "airr2.csv", """\
            sequence_id,sequence_aa,v_call,j_call
            s1,EVQLVESGG,IGHV3-23*01,IGHJ4*02
        """)
        results = abverse.reverse_translate(p, n_processes=1)
        assert results[0]["v_call"] == "IGHV3-23*01"
        assert results[0]["j_call"] == "IGHJ4*02"

    def test_mixed_annotated_and_unannotated(self, tmp_path):
        """CSV with some rows having gene calls, others not: MMseqs2 runs only for unannotated."""
        import abverse
        from abverse._reconstruct import _translate_nt
        from unittest.mock import patch, MagicMock
        import polars as pl

        p = _write(tmp_path, "mixed.csv", """\
            sequence_id,sequence_aa,v_call,j_call
            s1,EVQLVESGG,IGHV3-23*01,IGHJ4*02
            s2,QVQLVQ,,
        """)

        # s2 has null gene calls → MMseqs2 should be called for it
        with patch("abverse._pipeline.search_v_germline") as mock_v, \
             patch("abverse._pipeline.search_j_germline") as mock_j, \
             patch("abverse._pipeline.merge_vj_results") as mock_merge:
            mock_v.return_value = pl.DataFrame(schema={
                "query": pl.Utf8, "v_target": pl.Utf8, "v_nident": pl.Int64,
                "v_qstart": pl.Int64, "v_qend": pl.Int64, "v_tstart": pl.Int64,
                "v_tend": pl.Int64, "v_qseq": pl.Utf8, "v_tseq": pl.Utf8,
            })
            mock_j.return_value = pl.DataFrame(schema={
                "query": pl.Utf8, "j_target": pl.Utf8, "j_nident": pl.Int64,
                "j_qstart": pl.Int64, "j_qend": pl.Int64, "j_tstart": pl.Int64,
                "j_tend": pl.Int64, "j_qseq": pl.Utf8, "j_tseq": pl.Utf8,
            })
            mock_merge.return_value = pl.DataFrame({
                "query": ["s2"],
                "v_target": [None], "v_qstart": [None], "v_qend": [None], "v_tstart": [None],
                "j_target": [None], "j_qstart": [None], "j_qend": [None], "j_tstart": [None],
            })

            results = abverse.reverse_translate(p, n_processes=1)
            # V search was called (for s2 only)
            mock_v.assert_called_once()
            # Merge was called
            mock_merge.assert_called_once()

        assert len(results) == 2
        ids = {r.id for r in results}
        assert "s1" in ids and "s2" in ids

    def test_pairplex_csv_produces_two_chains(self, tmp_path):
        import abverse
        from abverse._reconstruct import _translate_nt

        header = "name,sequence_id:0,sequence_aa:0,v_call:0,j_call:0,locus:0,sequence_id:1,sequence_aa:1,v_call:1,j_call:1,locus:1"
        row = "pair1,pair1_H,EVQLVESGG,IGHV3-23*01,IGHJ4*02,IGH,pair1_L,QSVLTQ,IGLV2-14*01,IGLJ2*01,IGL"
        p = _write(tmp_path, "pairplex.csv", f"{header}\n{row}")

        results = abverse.reverse_translate(p, n_processes=1)
        assert len(results) == 2
        for r in results:
            assert _translate_nt(r.sequence) in ("EVQLVESGG", "QSVLTQ")

    def test_input_order_preserved_csv(self, tmp_path):
        import abverse
        from abverse._reconstruct import _translate_nt

        seqs = ["EVQLVES", "QVQLVQ", "MKPWST", "WGQGTL"]
        lines = "\n".join(f"s{i},{seq}" for i, seq in enumerate(seqs))
        p = _write(tmp_path, "order.csv", f"id,seq\n{lines}")
        results = abverse.reverse_translate(p, n_processes=1)
        assert len(results) == 4
        for i, r in enumerate(results):
            assert _translate_nt(r.sequence) == seqs[i]

    def test_empty_csv_returns_empty_list(self, tmp_path):
        import abverse
        p = _write(tmp_path, "empty.csv", "id,sequence\n")
        results = abverse.reverse_translate(p, n_processes=1)
        assert results == []
