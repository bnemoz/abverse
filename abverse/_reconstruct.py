# Per-sequence NT reconstruction from AA + germline alignment coordinates.
#
# All functions are pure (no I/O) and picklable for ProcessPoolExecutor.
#
# Coordinate convention (matching MMseqs2 easy-search, 0-indexed inclusive):
#   v_qstart, v_qend  — query AA positions (0-indexed, inclusive)
#   v_tstart, v_tend  — V germline AA positions (0-indexed, inclusive)
#   j_qstart, j_qend  — query AA positions RELATIVE TO POST-V SUBSEQUENCE (0-indexed inclusive)
#   j_tstart, j_tend  — J germline AA positions (0-indexed, inclusive)
#
# The reconstruction converts each AA position to an NT codon using:
#   - V/J region: optimal_codon(target_aa, germline_codon)
#   - 5'/CDR3/3' overhangs: fallback_codon(target_aa)

from __future__ import annotations

from typing import Optional

import abutils

from ._codons import fallback_codon, optimal_codon

__all__ = ["reconstruct_sequence", "_reconstruct_batch"]


def _get_germline_codon(
    nt_seq: str,
    frame: int,
    aa_pos: int,  # 0-indexed AA position in the germline alignment
    tstart: int,  # 0-indexed germline AA start of alignment
) -> Optional[str]:
    """Extract the germline codon at germline AA position (tstart + aa_pos).

    Returns None if the codon is truncated (edge of germline sequence).
    """
    germ_aa_idx = tstart + aa_pos
    nt_offset = (frame - 1) + germ_aa_idx * 3
    codon = nt_seq[nt_offset : nt_offset + 3]
    if len(codon) != 3:
        return None
    return codon


def reconstruct_sequence(
    seq_id: str,
    aa_seq: str,
    # V hit (may be None if no V assigned)
    v_call: Optional[str],
    v_qstart: Optional[int],
    v_qend: Optional[int],
    v_tstart: Optional[int],
    # J hit (may be None if no J assigned); coords relative to post-V region
    j_call: Optional[str],
    j_qstart: Optional[int],
    j_qend: Optional[int],
    j_tstart: Optional[int],
    # Germline sequence dicts and frame maps
    v_nt_seqs: dict[str, str],
    j_nt_seqs: dict[str, str],
    v_frame_map: dict[str, int],
    j_frame_map: dict[str, int],
) -> abutils.Sequence:
    """Reconstruct NT sequence for a single antibody AA sequence.

    Returns an abutils.Sequence with the reconstructed NT sequence and
    annotations: v_call, j_call, reconstruction_method.

    Raises ValueError if the input AA sequence contains a stop codon ('*').
    """
    if "*" in aa_seq:
        pos = aa_seq.index("*")
        raise ValueError(
            f"Stop codon at AA position {pos} in sequence '{seq_id}'"
        )

    n_aa = len(aa_seq)
    nt_out = ["NNN"] * n_aa  # initialise; will be overwritten position by position

    has_v = v_call is not None and v_qstart is not None and v_qend is not None and v_tstart is not None
    has_j = j_call is not None and j_qstart is not None and j_qend is not None and j_tstart is not None

    # ── Absolute J coordinates ───────────────────────────────────────────────
    # j_qstart is relative to aa_seq[v_qend+1:]; convert to absolute.
    j_abs_start: Optional[int] = None
    j_abs_end: Optional[int] = None
    if has_v and has_j:
        j_abs_start = (v_qend + 1) + j_qstart
        j_abs_end = (v_qend + 1) + j_qend
    elif has_j:
        j_abs_start = j_qstart
        j_abs_end = j_qend

    # ── Guard: V/J overlap ───────────────────────────────────────────────────
    # V takes priority; J starts after V end.
    if has_v and has_j and j_abs_start is not None and j_abs_start <= v_qend:
        j_abs_start = v_qend + 1
        if j_abs_start > j_abs_end:
            has_j = False

    # ── Build codon per position ─────────────────────────────────────────────
    v_nt = v_nt_seqs.get(v_call) if has_v else None
    v_frame = v_frame_map.get(v_call, 1) if has_v else 1

    j_nt = j_nt_seqs.get(j_call) if has_j else None
    j_frame = j_frame_map.get(j_call, 1) if has_j else 1

    method_parts: list[str] = []

    for i, aa in enumerate(aa_seq):
        # ── V region ────────────────────────────────────────────────────────
        if has_v and v_qstart <= i <= v_qend and v_nt is not None:
            aa_offset = i - v_qstart  # position within the aligned region
            germ_codon = _get_germline_codon(v_nt, v_frame, aa_offset, v_tstart)
            if germ_codon is not None:
                nt_out[i] = optimal_codon(aa, germ_codon)
            else:
                nt_out[i] = fallback_codon(aa)
            continue

        # ── J region ────────────────────────────────────────────────────────
        if has_j and j_abs_start is not None and j_abs_start <= i <= j_abs_end and j_nt is not None:
            aa_offset = i - j_abs_start  # position within the aligned J region
            germ_codon = _get_germline_codon(j_nt, j_frame, aa_offset, j_tstart)
            if germ_codon is not None:
                nt_out[i] = optimal_codon(aa, germ_codon)
            else:
                nt_out[i] = fallback_codon(aa)
            continue

        # ── 5' overhang / CDR3 / 3' overhang ────────────────────────────────
        nt_out[i] = fallback_codon(aa)

    nt_sequence = "".join(nt_out)

    # Determine reconstruction method annotation
    if has_v and has_j:
        method = "germline_vj"
    elif has_v:
        method = "germline_v_only"
    elif has_j:
        method = "germline_j_only"
    else:
        method = "codon_frequency"

    # Validate: translated output must match input AA
    translated = _translate_nt(nt_sequence)
    if translated != aa_seq:
        raise AssertionError(
            f"Reconstruction validation failed for '{seq_id}': "
            f"translate(output) != input_aa\n"
            f"  input : {aa_seq}\n"
            f"  output: {translated}"
        )

    result = abutils.Sequence(nt_sequence, id=seq_id)
    result["v_call"] = v_call
    result["j_call"] = j_call
    result["reconstruction_method"] = method
    return result


def _translate_nt(nt: str) -> str:
    """Translate a NT string in frame 1. Used only for validation."""
    from ._codons import CODON_TABLE
    aa_parts: list[str] = []
    for i in range(0, len(nt) - 2, 3):
        codon = nt[i : i + 3].upper()
        aa_parts.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa_parts)


def _reconstruct_batch(
    records: list[dict],
    v_nt_seqs: dict[str, str],
    j_nt_seqs: dict[str, str],
    v_frame_map: dict[str, int],
    j_frame_map: dict[str, int],
) -> list[abutils.Sequence | Exception]:
    """Process a batch of records; returns results or exceptions (one per record).

    Each record dict must have the keys expected by reconstruct_sequence.
    This function is top-level so it is picklable.
    """
    results: list[abutils.Sequence | Exception] = []
    for rec in records:
        try:
            seq = reconstruct_sequence(
                seq_id=rec["seq_id"],
                aa_seq=rec["aa_seq"],
                v_call=rec.get("v_call"),
                v_qstart=rec.get("v_qstart"),
                v_qend=rec.get("v_qend"),
                v_tstart=rec.get("v_tstart"),
                j_call=rec.get("j_call"),
                j_qstart=rec.get("j_qstart"),
                j_qend=rec.get("j_qend"),
                j_tstart=rec.get("j_tstart"),
                v_nt_seqs=v_nt_seqs,
                j_nt_seqs=j_nt_seqs,
                v_frame_map=v_frame_map,
                j_frame_map=j_frame_map,
            )
            results.append(seq)
        except Exception as exc:
            results.append(exc)
    return results
