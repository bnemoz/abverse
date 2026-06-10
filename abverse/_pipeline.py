# Orchestration pipeline: input normalization → germline DB → MMseqs2 V/J search
# → parallel per-sequence reconstruction → list[abutils.Sequence].

from __future__ import annotations

import os
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Iterable, Optional, Union

import abutils

from ._csv import parse_csv
from ._errors import ReverseTranslationError
from ._germline_db import build_germline_aa_db
from ._reconstruct import _reconstruct_batch
from ._search import (
    build_j_query_fasta,
    merge_vj_results,
    search_j_germline,
    search_v_germline,
)

__all__ = ["reverse_translate"]

SequenceInput = Union[str, abutils.Sequence, Iterable[abutils.Sequence]]

_CSV_EXTENSIONS = (".csv", ".tsv")


# ── Input normalisation ───────────────────────────────────────────────────────

def _normalise_input(
    sequences: SequenceInput,
) -> tuple[dict[str, str], dict[str, dict]]:
    """Return *(aa_seqs, pre_annotations)* from any supported input format.

    pre_annotations maps seq_id → dict with optional keys:
      v_call, j_call, v_qstart, v_qend, v_tstart, j_qstart, j_qend, j_tstart
    It is empty when the input carries no pre-existing annotation.
    """
    if isinstance(sequences, str):
        path = sequences
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Input file not found: {path}")
        if path.lower().endswith(_CSV_EXTENSIONS):
            return parse_csv(path)
        # FASTA
        seqs = abutils.io.read_fasta(path)
        return {s.id: str(s.sequence).upper() for s in seqs}, {}

    if isinstance(sequences, abutils.Sequence):
        sequences = [sequences]

    result: dict[str, str] = {}
    for i, s in enumerate(sequences):
        if isinstance(s, abutils.Sequence):
            sid = s.id or f"seq_{i}"
            seq = str(s.sequence).upper()
        elif isinstance(s, str):
            sid = f"seq_{i}"
            seq = s.upper()
        else:
            raise TypeError(f"Unsupported sequence type: {type(s)}")
        result[sid] = seq
    return result, {}


_STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _validate_residues(aa_seqs: dict[str, str]) -> list[dict]:
    """Return a list of ``invalid_residue`` failure dicts (empty if all valid).

    Only the 20 standard amino acids are accepted. Any other character
    (``X``, ``B``, ``Z``, ``U``, ``O``, ``*``, …) is reported. Positions in the
    detail message are 1-indexed for readability, grouped by offending residue.
    """
    failures: list[dict] = []
    for sid, seq in aa_seqs.items():
        bad: dict[str, list[int]] = {}
        for i, aa in enumerate(seq):
            if aa not in _STANDARD_AA:
                bad.setdefault(aa, []).append(i + 1)  # 1-indexed
        if not bad:
            continue
        parts: list[str] = []
        for residue, positions in bad.items():
            if len(positions) == 1:
                parts.append(f"'{residue}' at position {positions[0]}")
            else:
                joined = ", ".join(str(p) for p in positions)
                parts.append(f"'{residue}' at positions {joined}")
        failures.append({
            "seq_id": sid,
            "kind": "invalid_residue",
            "detail": "; ".join(parts),
        })
    return failures


def _write_aa_fasta(aa_seqs: dict[str, str], path: str) -> None:
    with open(path, "w") as fh:
        for sid, seq in aa_seqs.items():
            fh.write(f">{sid}\n{seq}\n")


# ── Record builders ───────────────────────────────────────────────────────────

def _df_to_records(merged_df, aa_seqs: dict[str, str]) -> list[dict]:
    records: list[dict] = []
    for row in merged_df.iter_rows(named=True):
        sid = row["query"]
        rec: dict = {
            "seq_id": sid,
            "aa_seq": aa_seqs[sid],
            "v_call": row.get("v_target"),
            "v_qstart": row.get("v_qstart"),
            "v_qend": row.get("v_qend"),
            "v_tstart": row.get("v_tstart"),
            "j_call": row.get("j_target"),
            "j_qstart": row.get("j_qstart"),
            "j_qend": row.get("j_qend"),
            "j_tstart": row.get("j_tstart"),
        }
        records.append(rec)
    return records


def _annotations_to_records(
    aa_seqs: dict[str, str],
    pre_annotations: dict[str, dict],
    seq_ids: set[str],
) -> list[dict]:
    """Build reconstruction records directly from pre-supplied annotations."""
    records: list[dict] = []
    for sid in seq_ids:
        ann = pre_annotations.get(sid, {})
        records.append({
            "seq_id": sid,
            "aa_seq": aa_seqs[sid],
            "v_call":   ann.get("v_call"),
            "v_qstart": ann.get("v_qstart"),
            "v_qend":   ann.get("v_qend"),
            "v_tstart": ann.get("v_tstart"),
            "j_call":   ann.get("j_call"),
            "j_qstart": ann.get("j_qstart"),
            "j_qend":   ann.get("j_qend"),
            "j_tstart": ann.get("j_tstart"),
        })
    return records


# ── Parallel dispatch ─────────────────────────────────────────────────────────

def _chunk_records(records: list[dict], chunksize: int) -> list[list[dict]]:
    return [records[i : i + chunksize] for i in range(0, len(records), chunksize)]


# ── Public API ────────────────────────────────────────────────────────────────

def reverse_translate(
    sequences: SequenceInput,
    species: str = "human",
    receptor: str = "bcr",
    n_processes: Optional[int] = None,
    threads: Optional[int] = None,
    chunksize: int = 500,
    force_rebuild_db: bool = False,
    output_fasta: Optional[str] = None,
    verbose: bool = False,
) -> list[abutils.Sequence]:
    """Reverse-translate antibody AA sequences to germline-faithful NT sequences.

    Parameters
    ----------
    sequences:
        Any of the following:

        * **FASTA file path** (``.fasta`` / ``.fa``): one sequence per record.
        * **CSV / TSV file path**: auto-detected format —

          - *Simple CSV*: first column = sequence ID, second column = AA sequence.
          - *AIRR TSV/CSV*: AIRR-standard columns (``sequence_id``,
            ``sequence_aa``, ``v_call``, ``j_call``, …).  If germline gene
            calls are present the MMseqs2 search is skipped for those
            sequences; if AIRR alignment coordinates are also present they
            are used directly for germline-informed codon selection.
          - *PairPlex TSV/CSV*: paired output with ``:0`` / ``:1`` column
            suffixes.  Each antibody row is split into two chain records.
            Pre-supplied gene calls and coordinates are handled identically
            to the AIRR path.

        * **abutils.Sequence** (single or iterable).
        * **Iterable of plain strings** (sequences only; IDs auto-assigned).
    species:
        Germline species.  Currently only 'human' is supported.
    receptor:
        Receptor type.  Currently only 'bcr' is supported.
    n_processes:
        Number of worker processes for parallel reconstruction.
        Defaults to os.cpu_count().
    threads:
        Threads passed to MMseqs2.  Defaults to MMseqs2's auto-detection.
    chunksize:
        Number of sequences per worker batch.
    force_rebuild_db:
        Force re-translation and re-build of germline AA databases.
    output_fasta:
        Optional path to write the reconstructed NT FASTA.
    verbose:
        Print progress messages.

    Returns
    -------
    list[abutils.Sequence]
        Reconstructed NT sequences in the same order as the input.

    Raises
    ------
    ReverseTranslationError
        If any input sequence contains a residue outside the 20 standard
        amino acids, or if reconstruction fails for any sequence. All failures
        are collected and reported together; the error's ``.failures`` attribute
        lists each offending sequence id, kind, and detail. Input is validated
        before any germline/MMseqs2 work, so bad input fails fast.
    """
    if n_processes is None:
        n_processes = os.cpu_count() or 1

    def _log(msg: str) -> None:
        if verbose:
            print(f"[abverse] {msg}")

    # 1. Normalise input
    _log("Normalising input sequences …")
    aa_seqs, pre_annotations = _normalise_input(sequences)
    if not aa_seqs:
        return []

    # Validate residues upfront across every input path; collect all failures
    # and raise one comprehensive error before doing any germline/MMseqs2 work.
    residue_failures = _validate_residues(aa_seqs)
    if residue_failures:
        raise ReverseTranslationError(residue_failures)

    input_order = list(aa_seqs.keys())
    _log(f"  {len(aa_seqs):,} sequences")

    # Separate sequences: those with pre-supplied gene calls skip MMseqs2.
    annotated_ids = {
        sid for sid, ann in pre_annotations.items()
        if ann.get("v_call") is not None or ann.get("j_call") is not None
    }
    unannotated_ids = {sid for sid in aa_seqs if sid not in annotated_ids}

    if annotated_ids:
        _log(
            f"  {len(annotated_ids):,} sequences have pre-supplied gene calls "
            f"(MMseqs2 search skipped for these)"
        )

    with tempfile.TemporaryDirectory(prefix="abverse_") as tmpdir:
        # 2. Build / load germline AA databases (always needed for NT codon lookup)
        _log("Loading germline AA databases …")
        db = build_germline_aa_db(species=species, receptor=receptor, force_rebuild=force_rebuild_db)

        records: list[dict] = []

        # 3a. Records from pre-annotations (no MMseqs2)
        if annotated_ids:
            records.extend(
                _annotations_to_records(aa_seqs, pre_annotations, annotated_ids)
            )

        # 3b. Records from MMseqs2 search for unannotated sequences
        if unannotated_ids:
            unannotated_seqs = {sid: aa_seqs[sid] for sid in unannotated_ids}

            # Write AA query FASTA
            aa_fasta = os.path.join(tmpdir, "query_aa.fasta")
            _write_aa_fasta(unannotated_seqs, aa_fasta)

            # V search
            _log("Running V germline search …")
            v_result_path = os.path.join(tmpdir, "v_results.tsv")
            v_df = search_v_germline(
                query_fasta=aa_fasta,
                v_db_path=db["v_db_path"],
                output_path=v_result_path,
                threads=threads,
            )
            _log(f"  {len(v_df):,} V assignments")

            # Build J query FASTA (post-V region)
            j_query_fasta = os.path.join(tmpdir, "j_query.fasta")
            build_j_query_fasta(
                v_results=v_df,
                input_aa_seqs=unannotated_seqs,
                output_path=j_query_fasta,
            )

            # J search
            _log("Running J germline search …")
            j_result_path = os.path.join(tmpdir, "j_results.tsv")
            j_df = search_j_germline(
                j_query_fasta=j_query_fasta,
                j_db_path=db["j_db_path"],
                output_path=j_result_path,
                threads=threads,
            )
            _log(f"  {len(j_df):,} J assignments")

            merged_df = merge_vj_results(v_df, j_df, unannotated_seqs)
            records.extend(_df_to_records(merged_df, unannotated_seqs))

        # 4. Dispatch reconstruction in parallel
        _log(f"Reconstructing NT sequences (workers={n_processes}) …")
        chunks = _chunk_records(records, chunksize)

        results_by_id: dict[str, abutils.Sequence] = {}
        # Reconstruction failures are collected and raised together. Post input
        # validation these should essentially never fire for legitimate input,
        # so a raise here surfaces a genuine bug immediately rather than emitting
        # junk that fails downstream.
        recon_failures: list[dict] = []

        if n_processes == 1 or len(records) <= chunksize:
            # Single-process path (avoids spawn overhead for small inputs)
            batch_results = _reconstruct_batch(
                records,
                db["v_nt_seqs"],
                db["j_nt_seqs"],
                db["v_frame_map"],
                db["j_frame_map"],
            )
            for rec, res in zip(records, batch_results):
                if isinstance(res, Exception):
                    recon_failures.append({
                        "seq_id": rec["seq_id"],
                        "kind": "reconstruction_error",
                        "detail": str(res),
                    })
                else:
                    results_by_id[res.id] = res
        else:
            with ProcessPoolExecutor(
                max_workers=n_processes,
                mp_context=__import__("multiprocessing").get_context("spawn"),
            ) as executor:
                futures = {
                    executor.submit(
                        _reconstruct_batch,
                        chunk,
                        db["v_nt_seqs"],
                        db["j_nt_seqs"],
                        db["v_frame_map"],
                        db["j_frame_map"],
                    ): chunk
                    for chunk in chunks
                }
                for future in as_completed(futures):
                    chunk = futures[future]
                    try:
                        batch_results = future.result()
                    except Exception as exc:
                        # Entire chunk failed — record every record in it
                        for rec in chunk:
                            recon_failures.append({
                                "seq_id": rec["seq_id"],
                                "kind": "reconstruction_error",
                                "detail": str(exc),
                            })
                        continue
                    for rec, res in zip(chunk, batch_results):
                        if isinstance(res, Exception):
                            recon_failures.append({
                                "seq_id": rec["seq_id"],
                                "kind": "reconstruction_error",
                                "detail": str(res),
                            })
                        else:
                            results_by_id[res.id] = res

        if recon_failures:
            raise ReverseTranslationError(recon_failures)

    # 5. Return in input order
    _log("Done.")
    ordered = [results_by_id[sid] for sid in input_order if sid in results_by_id]

    # 6. Optional FASTA write
    if output_fasta is not None:
        abutils.io.to_fasta(ordered, output_fasta)

    return ordered
