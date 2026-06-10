"""Tests for _reconstruct.py — per-sequence NT reconstruction logic."""

import pytest
import abutils

from abverse._codons import CODON_TABLE, fallback_codon
from abverse._reconstruct import _translate_nt, reconstruct_sequence


def _build_simple_db():
    """Return minimal v_nt_seqs / j_nt_seqs for test purposes."""
    # A V gene: 3 AA = 9 nt, frame 1
    # Encodes QVQ
    v_nt = "CAGGTGCAG"  # Q V Q
    # A J gene: 3 AA = 9 nt, frame 1
    # Encodes WGQ
    j_nt = "TGGGGACAG"  # W G Q
    return (
        {"IGHV1-1*01": v_nt},
        {"IGHJ1-1*01": j_nt},
        {"IGHV1-1*01": 1},
        {"IGHJ1-1*01": 1},
    )


class TestTranslateNt:
    def test_basic(self):
        assert _translate_nt("ATGAAACCC") == "MKP"

    def test_stop_codon(self):
        assert _translate_nt("ATGTAA") == "M*"

    def test_truncated(self):
        # 1 spare nt — ignored
        assert _translate_nt("ATGA") == "M"


class TestReconstructSequence:
    def test_simple_v_and_j(self):
        v_nt_seqs, j_nt_seqs, v_frame_map, j_frame_map = _build_simple_db()
        # AA: Q V Q A R W G Q   (A R = CDR3, standard residues)
        # V covers pos 0..2 (QVQ), J covers pos 5..7 relative to post-V = abs 3..5 wait…
        # post-V starts at v_qend+1 = 3; j_qstart=0 → abs j_start = 3
        aa_seq = "QVQARWGQ"
        result = reconstruct_sequence(
            seq_id="test",
            aa_seq=aa_seq,
            v_call="IGHV1-1*01",
            v_qstart=0,
            v_qend=2,
            v_tstart=0,
            j_call="IGHJ1-1*01",
            j_qstart=2,  # relative to post-V (aa_seq[3:]) → abs 5
            j_qend=4,
            j_tstart=0,
            v_nt_seqs=v_nt_seqs,
            j_nt_seqs=j_nt_seqs,
            v_frame_map=v_frame_map,
            j_frame_map=j_frame_map,
        )
        assert isinstance(result, abutils.Sequence)
        assert len(result.sequence) == len(aa_seq) * 3
        # Translation must match input
        assert _translate_nt(result.sequence) == aa_seq

    def test_translation_matches_input(self):
        v_nt_seqs, j_nt_seqs, v_frame_map, j_frame_map = _build_simple_db()
        aa_seq = "QVQMKP"
        result = reconstruct_sequence(
            seq_id="t2",
            aa_seq=aa_seq,
            v_call="IGHV1-1*01",
            v_qstart=0,
            v_qend=2,
            v_tstart=0,
            j_call=None,
            j_qstart=None,
            j_qend=None,
            j_tstart=None,
            v_nt_seqs=v_nt_seqs,
            j_nt_seqs=j_nt_seqs,
            v_frame_map=v_frame_map,
            j_frame_map=j_frame_map,
        )
        assert _translate_nt(result.sequence) == aa_seq

    def test_no_v_no_j_uses_fallback(self):
        v_nt_seqs, j_nt_seqs, v_frame_map, j_frame_map = _build_simple_db()
        aa_seq = "MKPWST"
        result = reconstruct_sequence(
            seq_id="t3",
            aa_seq=aa_seq,
            v_call=None,
            v_qstart=None,
            v_qend=None,
            v_tstart=None,
            j_call=None,
            j_qstart=None,
            j_qend=None,
            j_tstart=None,
            v_nt_seqs=v_nt_seqs,
            j_nt_seqs=j_nt_seqs,
            v_frame_map=v_frame_map,
            j_frame_map=j_frame_map,
        )
        assert _translate_nt(result.sequence) == aa_seq
        assert result["reconstruction_method"] == "codon_frequency"

    def test_stop_codon_raises(self):
        v_nt_seqs, j_nt_seqs, v_frame_map, j_frame_map = _build_simple_db()
        with pytest.raises(ValueError, match="Stop codon"):
            reconstruct_sequence(
                seq_id="stop",
                aa_seq="MK*P",
                v_call=None, v_qstart=None, v_qend=None, v_tstart=None,
                j_call=None, j_qstart=None, j_qend=None, j_tstart=None,
                v_nt_seqs=v_nt_seqs,
                j_nt_seqs=j_nt_seqs,
                v_frame_map=v_frame_map,
                j_frame_map=j_frame_map,
            )

    def test_nonstandard_aa_raises(self):
        # Reconstruction assumes input is already validated; a non-standard
        # residue reaching this layer is a genuine error and must raise, not
        # silently pass through as NNN.
        v_nt_seqs, j_nt_seqs, v_frame_map, j_frame_map = _build_simple_db()
        with pytest.raises((KeyError, ValueError)):
            reconstruct_sequence(
                seq_id="nonstandard",
                aa_seq="MXK",
                v_call=None, v_qstart=None, v_qend=None, v_tstart=None,
                j_call=None, j_qstart=None, j_qend=None, j_tstart=None,
                v_nt_seqs=v_nt_seqs,
                j_nt_seqs=j_nt_seqs,
                v_frame_map=v_frame_map,
                j_frame_map=j_frame_map,
            )

    def test_vj_overlap_v_priority(self):
        """When J abs start overlaps with V end, V takes priority and J begins after V."""
        v_nt_seqs, j_nt_seqs, v_frame_map, j_frame_map = _build_simple_db()
        aa_seq = "QVQWGQ"
        # J j_qstart=0 relative to post-V (aa[3:]) → abs start = 3; v_qend = 2 → no overlap
        # Force overlap: j_qstart = -2 (abs = 3 + (-2) = 1) — overlap with V (v_qend=2)
        # The code converts j_abs_start = (v_qend+1) + j_qstart = 3 + (-2) = 1 → overlap → fix to 3
        result = reconstruct_sequence(
            seq_id="overlap",
            aa_seq=aa_seq,
            v_call="IGHV1-1*01",
            v_qstart=0,
            v_qend=2,
            v_tstart=0,
            j_call="IGHJ1-1*01",
            j_qstart=-2,  # forces overlap
            j_qend=1,
            j_tstart=0,
            v_nt_seqs=v_nt_seqs,
            j_nt_seqs=j_nt_seqs,
            v_frame_map=v_frame_map,
            j_frame_map=j_frame_map,
        )
        assert _translate_nt(result.sequence) == aa_seq

    def test_annotations_set(self):
        v_nt_seqs, j_nt_seqs, v_frame_map, j_frame_map = _build_simple_db()
        result = reconstruct_sequence(
            seq_id="annot",
            aa_seq="QVQWGQ",
            v_call="IGHV1-1*01",
            v_qstart=0, v_qend=2, v_tstart=0,
            j_call="IGHJ1-1*01",
            j_qstart=0, j_qend=2, j_tstart=0,
            v_nt_seqs=v_nt_seqs,
            j_nt_seqs=j_nt_seqs,
            v_frame_map=v_frame_map,
            j_frame_map=j_frame_map,
        )
        assert result["v_call"] == "IGHV1-1*01"
        assert result["j_call"] == "IGHJ1-1*01"
        assert result["reconstruction_method"] == "germline_vj"

    def test_germline_codon_at_edge_uses_fallback(self):
        """Germline codon truncated at edge of germline sequence → fallback."""
        v_nt_seqs, j_nt_seqs, v_frame_map, j_frame_map = _build_simple_db()
        aa_seq = "QVQ"
        # v_tstart=1 so position 2 → germline AA idx 3 → nt offset 9 = beyond 9-nt germline
        result = reconstruct_sequence(
            seq_id="edge",
            aa_seq=aa_seq,
            v_call="IGHV1-1*01",
            v_qstart=0, v_qend=2, v_tstart=1,
            j_call=None, j_qstart=None, j_qend=None, j_tstart=None,
            v_nt_seqs=v_nt_seqs,
            j_nt_seqs=j_nt_seqs,
            v_frame_map=v_frame_map,
            j_frame_map=j_frame_map,
        )
        assert _translate_nt(result.sequence) == aa_seq
