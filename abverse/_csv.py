# CSV input parsing for abverse.
#
# Supports three formats, auto-detected from column names:
#   simple   — col 0 = seq_id, col 1 = aa_sequence  (no pre-annotation)
#   airr     — AIRR-standard columns (sequence_id, sequence_aa, v_call, j_call, …)
#   pairplex — PairPlex paired output: one row per antibody, columns suffixed :0 / :1
#
# Returns (aa_seqs, pre_annotations):
#   aa_seqs          : {seq_id: aa_sequence}
#   pre_annotations  : {seq_id: {v_call, j_call, v_qstart, v_qend, v_tstart,
#                                j_qstart, j_qend, j_tstart}}
#
# pre_annotations is empty ({}) for simple CSVs.
# For AIRR / PairPlex, annotation keys present only when the corresponding
# columns exist in the file; coordinate keys present only when all six AIRR
# coordinate columns are populated for that row.
#
# AIRR coordinate conversion (1-indexed NT → 0-indexed AA):
#   aa_pos = (nt_pos_1indexed - 1) // 3
# Assumption: sequence_aa begins at the first codon of the V gene
# (v_sequence_start = 1 for most productive BCR / TCR sequences).
# j_qstart / j_qend are stored relative to the post-V sub-sequence,
# matching the convention in _reconstruct.py.

from __future__ import annotations

import os
from typing import Literal, Optional

import polars as pl

__all__ = ["parse_csv", "detect_csv_format"]

# AIRR columns required to attempt coordinate extraction
_AIRR_COORD_COLS = {
    "v_sequence_start",
    "v_sequence_end",
    "v_germline_start",
    "j_sequence_start",
    "j_sequence_end",
    "j_germline_start",
}


# ── Format detection ──────────────────────────────────────────────────────────

def detect_csv_format(
    columns: list[str],
) -> Literal["simple", "airr", "pairplex"]:
    """Infer CSV format from column names."""
    col_set = set(columns)
    if any(c.endswith(":0") or c.endswith(":1") for c in col_set):
        return "pairplex"
    if "sequence_id" in col_set and (
        "sequence_aa" in col_set or "sequence" in col_set
    ):
        return "airr"
    return "simple"


# ── Public entry point ────────────────────────────────────────────────────────

def parse_csv(path: str) -> tuple[dict[str, str], dict[str, dict]]:
    """Parse a CSV or TSV file and return *(aa_seqs, pre_annotations)*.

    Parameters
    ----------
    path:
        Path to the file.  Tab-separated if the extension is ``.tsv``,
        comma-separated otherwise.

    Returns
    -------
    aa_seqs : dict[str, str]
        ``{seq_id: aa_sequence}``
    pre_annotations : dict[str, dict]
        ``{seq_id: {v_call, j_call, [v_qstart, v_qend, v_tstart,
        j_qstart, j_qend, j_tstart]}}``
        Empty when no annotations are available (simple CSV).
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"CSV file not found: {path}")

    sep = "\t" if path.endswith(".tsv") else ","
    df = pl.read_csv(
        path,
        separator=sep,
        null_values=["", "NA", "N/A", "None", "none", "nan", "NaN"],
        infer_schema_length=10_000,
    )

    fmt = detect_csv_format(df.columns)

    if fmt == "pairplex":
        return _parse_pairplex(df)
    if fmt == "airr":
        return _parse_airr(df)
    return _parse_simple(df)


# ── Simple CSV ────────────────────────────────────────────────────────────────

def _parse_simple(df: pl.DataFrame) -> tuple[dict[str, str], dict[str, dict]]:
    """Two-column CSV: first col = seq_id, second col = AA sequence."""
    if df.width < 2:
        raise ValueError(
            "Simple CSV must have at least 2 columns (seq_id, aa_sequence)"
        )
    id_col, seq_col = df.columns[0], df.columns[1]
    aa_seqs: dict[str, str] = {}
    for row in df.iter_rows(named=True):
        sid = row[id_col]
        seq = row[seq_col]
        if sid is None or seq is None:
            continue
        sid = str(sid)
        seq = str(seq).upper().strip()
        if seq:
            aa_seqs[sid] = seq
    return aa_seqs, {}


# ── AIRR ──────────────────────────────────────────────────────────────────────

def _parse_airr(df: pl.DataFrame) -> tuple[dict[str, str], dict[str, dict]]:
    """Parse AIRR-standard CSV/TSV (unpaired, one row per chain)."""
    cols = set(df.columns)

    if "sequence_aa" not in cols:
        raise ValueError(
            "AIRR CSV must contain a 'sequence_aa' column.  "
            "If only a nucleotide 'sequence' column is present, "
            "translate it before calling abverse."
        )

    has_gene_calls = {"v_call", "j_call"}.issubset(cols)
    has_coords = _AIRR_COORD_COLS.issubset(cols)

    aa_seqs: dict[str, str] = {}
    pre_annotations: dict[str, dict] = {}

    for row in df.iter_rows(named=True):
        sid = row.get("sequence_id")
        seq = row.get("sequence_aa")
        if sid is None or seq is None:
            continue
        sid = str(sid)
        seq = _clean_aa(str(seq))
        if not seq:
            continue
        aa_seqs[sid] = seq

        if has_gene_calls:
            annot: dict = {
                "v_call": row.get("v_call"),
                "j_call": row.get("j_call"),
            }
            if has_coords:
                coords = _extract_airr_coords(row)
                annot.update(coords)
            pre_annotations[sid] = annot

    return aa_seqs, pre_annotations


# ── PairPlex ──────────────────────────────────────────────────────────────────

def _parse_pairplex(df: pl.DataFrame) -> tuple[dict[str, str], dict[str, dict]]:
    """Parse PairPlex paired output (columns suffixed ``:0`` / ``:1``).

    Each row represents one antibody; ``:0`` is the first chain (typically
    heavy, ``locus:0 == IGH``) and ``:1`` is the second chain (light).
    The pair identifier is taken from the ``name`` column when present,
    otherwise from the first column.
    """
    cols = set(df.columns)
    pair_id_col = "name" if "name" in cols else df.columns[0]

    aa_seqs: dict[str, str] = {}
    pre_annotations: dict[str, dict] = {}

    for suffix in ("0", "1"):
        seq_aa_col = f"sequence_aa:{suffix}"
        if seq_aa_col not in cols:
            continue

        seq_id_col = f"sequence_id:{suffix}" if f"sequence_id:{suffix}" in cols else None
        v_call_col = f"v_call:{suffix}" if f"v_call:{suffix}" in cols else None
        j_call_col = f"j_call:{suffix}" if f"j_call:{suffix}" in cols else None
        locus_col  = f"locus:{suffix}"  if f"locus:{suffix}"  in cols else None

        coord_cols_present = all(
            f"{base}:{suffix}" in cols
            for base in (
                "v_sequence_start", "v_sequence_end", "v_germline_start",
                "j_sequence_start", "j_sequence_end", "j_germline_start",
            )
        )

        for row in df.iter_rows(named=True):
            seq = row.get(seq_aa_col)
            if seq is None:
                continue
            seq = _clean_aa(str(seq))
            if not seq:
                continue

            # Build seq_id: prefer explicit sequence_id:N, else name + chain label
            if seq_id_col and row.get(seq_id_col):
                sid = str(row[seq_id_col])
            else:
                pair_id = str(row[pair_id_col])
                label = _chain_label(locus_col, row, suffix)
                sid = f"{pair_id}_{label}"

            aa_seqs[sid] = seq

            annot: dict = {}
            if v_call_col:
                annot["v_call"] = row.get(v_call_col)
            if j_call_col:
                annot["j_call"] = row.get(j_call_col)

            if coord_cols_present:
                pseudo_row = {
                    "v_sequence_start": row.get(f"v_sequence_start:{suffix}"),
                    "v_sequence_end":   row.get(f"v_sequence_end:{suffix}"),
                    "v_germline_start": row.get(f"v_germline_start:{suffix}"),
                    "j_sequence_start": row.get(f"j_sequence_start:{suffix}"),
                    "j_sequence_end":   row.get(f"j_sequence_end:{suffix}"),
                    "j_germline_start": row.get(f"j_germline_start:{suffix}"),
                }
                annot.update(_extract_airr_coords(pseudo_row))

            if annot:
                pre_annotations[sid] = annot

    return aa_seqs, pre_annotations


# ── Coordinate extraction ─────────────────────────────────────────────────────

def _extract_airr_coords(row: dict) -> dict:
    """Convert AIRR 1-indexed NT coordinates to 0-indexed AA coordinates.

    Returns a (possibly partial) dict with keys:
      v_qstart, v_qend, v_tstart, j_qstart, j_qend, j_tstart

    Missing or None source values silently omit the corresponding keys.
    j_qstart / j_qend are stored relative to the post-V subsequence so that
    they match the convention expected by reconstruct_sequence().
    """
    result: dict = {}

    vs  = _to_int(row.get("v_sequence_start"))
    ve  = _to_int(row.get("v_sequence_end"))
    vgs = _to_int(row.get("v_germline_start"))

    v_qend: Optional[int] = None
    if vs is not None and ve is not None and vgs is not None:
        result["v_qstart"] = (vs - 1) // 3
        result["v_qend"]   = (ve - 1) // 3
        result["v_tstart"] = (vgs - 1) // 3
        v_qend = result["v_qend"]

    js  = _to_int(row.get("j_sequence_start"))
    je  = _to_int(row.get("j_sequence_end"))
    jgs = _to_int(row.get("j_germline_start"))

    if js is not None and je is not None and jgs is not None:
        j_abs_start = (js - 1) // 3
        j_abs_end   = (je - 1) // 3
        result["j_tstart"] = (jgs - 1) // 3

        if v_qend is not None:
            result["j_qstart"] = j_abs_start - (v_qend + 1)
            result["j_qend"]   = j_abs_end   - (v_qend + 1)
        else:
            result["j_qstart"] = j_abs_start
            result["j_qend"]   = j_abs_end

    return result


# ── Helpers ───────────────────────────────────────────────────────────────────

def _clean_aa(seq: str) -> str:
    """Strip whitespace and IMGT gap characters; upper-case."""
    return seq.replace(".", "").replace("-", "").upper().strip()


def _to_int(value) -> Optional[int]:
    """Return int or None (handles None, NaN, and string representations)."""
    if value is None:
        return None
    try:
        return int(value)
    except (ValueError, TypeError):
        return None


def _chain_label(locus_col: Optional[str], row: dict, suffix: str) -> str:
    """Return 'H' for IGH locus, 'L' for IGK/IGL, else 'chain{suffix}'."""
    if locus_col:
        locus = row.get(locus_col)
        if locus == "IGH":
            return "H"
        if locus in ("IGK", "IGL"):
            return "L"
    return f"chain{suffix}"
