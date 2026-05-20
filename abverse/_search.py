# MMseqs2 protein-protein search wrappers for V and J germline assignment.
#
# V search: high sensitivity, protein-protein (search_type=1).
# J search: looser E-value threshold because J genes are short (~15 AA).
# Results parsed with Polars; best hit per query kept by nident (descending).
#
# Coordinate convention (MMseqs2 easy-search, 0-indexed inclusive on both ends):
#   qstart, qend  — query positions
#   tstart, tend  — target (germline) positions
# v_qend is the last aligned query AA (inclusive).
# j_qstart/j_qend are relative to the post-V subsequence passed to j_search.

from __future__ import annotations

import os
from typing import Optional

import abutils
import abutils.tl as tl
import polars as pl

__all__ = [
    "search_v_germline",
    "build_j_query_fasta",
    "search_j_germline",
    "merge_vj_results",
]

_FORMAT = "query,target,nident,qstart,qend,tstart,tend,qseq,tseq"


def search_v_germline(
    query_fasta: str,
    v_db_path: str,
    output_path: str,
    threads: Optional[int] = None,
) -> pl.DataFrame:
    """Run protein-protein V search and return a DataFrame of best hits."""
    tl.mmseqs_search(
        query=query_fasta,
        target=v_db_path,
        output_path=output_path,
        search_type=1,
        max_seqs=25,
        sensitivity=7.5,
        format_mode=4,
        format_output=_FORMAT,
        additional_cli_args="--min-aln-len 10 --alignment-mode 3",
        threads=threads,
        verbosity=0,
    )
    return _parse_results(output_path, prefix="v")


def build_j_query_fasta(
    v_results: pl.DataFrame,
    input_aa_seqs: dict[str, str],
    output_path: str,
) -> None:
    """Write a FASTA with aa_seq[v_qend+1:] per sequence (post-V region for J search).

    Sequences with no V hit get the full AA sequence.
    Sequences whose V alignment reaches the end are skipped (no J region possible).
    """
    v_map = {}
    for row in v_results.iter_rows(named=True):
        v_map[row["query"]] = row["v_qend"]  # 0-indexed inclusive

    with open(output_path, "w") as fh:
        for seq_id, aa_seq in input_aa_seqs.items():
            v_qend = v_map.get(seq_id, -1)
            post_v = aa_seq[v_qend + 1 :]
            if len(post_v) < 3:
                continue
            fh.write(f">{seq_id}\n{post_v}\n")


def search_j_germline(
    j_query_fasta: str,
    j_db_path: str,
    output_path: str,
    threads: Optional[int] = None,
) -> pl.DataFrame:
    """Run protein-protein J search (loose E-value) and return best hits."""
    if not os.path.isfile(j_query_fasta) or os.path.getsize(j_query_fasta) == 0:
        return pl.DataFrame(
            schema={
                "query": pl.Utf8,
                "j_target": pl.Utf8,
                "j_nident": pl.Int64,
                "j_qstart": pl.Int64,
                "j_qend": pl.Int64,
                "j_tstart": pl.Int64,
                "j_tend": pl.Int64,
                "j_qseq": pl.Utf8,
                "j_tseq": pl.Utf8,
            }
        )
    tl.mmseqs_search(
        query=j_query_fasta,
        target=j_db_path,
        output_path=output_path,
        search_type=1,
        max_seqs=25,
        max_evalue=1000.0,
        format_mode=4,
        format_output=_FORMAT,
        additional_cli_args="--min-aln-len 5 --alignment-mode 3",
        threads=threads,
        verbosity=0,
    )
    return _parse_results(output_path, prefix="j")


def merge_vj_results(
    v_df: pl.DataFrame,
    j_df: pl.DataFrame,
    input_aa_seqs: dict[str, str],
    v_qend_map: Optional[dict[str, int]] = None,
) -> pl.DataFrame:
    """Merge V and J results into a single DataFrame indexed by sequence ID.

    j_qstart / j_qend are stored as-is (relative to post-V query).
    The reconstruction step converts them to absolute AA coordinates.

    Parameters
    ----------
    v_qend_map:
        Pre-computed {seq_id: v_qend} map.  If None, derived from v_df.
    """
    all_ids = list(input_aa_seqs.keys())
    base = pl.DataFrame({"query": all_ids})

    merged = base.join(v_df, on="query", how="left")

    if not j_df.is_empty():
        merged = merged.join(j_df, on="query", how="left")
    else:
        for col in ["j_target", "j_nident", "j_qstart", "j_qend", "j_tstart", "j_tend", "j_qseq", "j_tseq"]:
            merged = merged.with_columns(pl.lit(None).alias(col))

    return merged


# ── Internal helpers ──────────────────────────────────────────────────────────

def _parse_results(output_path: str, prefix: str) -> pl.DataFrame:
    """Parse MMseqs2 TSV output (format_mode=4 includes headers)."""
    if not os.path.isfile(output_path) or os.path.getsize(output_path) == 0:
        return pl.DataFrame(
            schema={
                "query": pl.Utf8,
                f"{prefix}_target": pl.Utf8,
                f"{prefix}_nident": pl.Int64,
                f"{prefix}_qstart": pl.Int64,
                f"{prefix}_qend": pl.Int64,
                f"{prefix}_tstart": pl.Int64,
                f"{prefix}_tend": pl.Int64,
                f"{prefix}_qseq": pl.Utf8,
                f"{prefix}_tseq": pl.Utf8,
            }
        )

    df = pl.read_csv(output_path, separator="\t", has_header=True)

    # Rename non-query columns to {prefix}_{col}
    rename_map = {c: f"{prefix}_{c}" if c != "query" else c for c in df.columns}
    df = df.rename(rename_map)

    # Keep best hit per query (highest nident)
    df = (
        df.sort(f"{prefix}_nident", descending=True, nulls_last=True)
        .unique(subset=["query"], keep="first")
    )
    return df
