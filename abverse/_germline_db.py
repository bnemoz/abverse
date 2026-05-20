# Build and cache amino-acid germline databases for MMseqs2 protein-protein search.
#
# Source NT germlines (from abstar) are translated to AA FASTA files, then converted
# to MMseqs2 protein databases. Results are cached under ~/.abverse/germline_dbs/ and
# invalidated via SHA-256 of the source NT FASTA.

from __future__ import annotations

import hashlib
import json
import os
import re
from typing import Optional

import abutils
import abutils.tl as tl

__all__ = ["build_germline_aa_db", "get_germline_aa_db_path"]

_CACHE_ROOT = os.path.expanduser("~/.abverse/germline_dbs")

# Conserved J-gene motifs used for frame detection.
# IGH J genes end with WG[^*]G; IGK/IGL end with FG[^*]G.
_J_MOTIF_RE = re.compile(r"[WF]G.G", re.IGNORECASE)


# ── Frame detection for J genes ──────────────────────────────────────────────

def _detect_j_frame(nt_seq: str) -> int:
    """Return reading frame (1, 2, or 3) for a J-gene NT sequence.

    Tries each frame and looks for the conserved WG.G / FG.G motif.
    Falls back to the frame that gives the longest stop-codon-free translation.
    """
    nt_seq = nt_seq.upper()
    for frame in (1, 2, 3):
        offset = frame - 1
        aa = _translate_nt(nt_seq[offset:])
        if _J_MOTIF_RE.search(aa):
            return frame
    # Fallback: longest stop-free translation
    best_frame = 1
    best_len = -1
    for frame in (1, 2, 3):
        offset = frame - 1
        aa = _translate_nt(nt_seq[offset:])
        clean = aa.split("*")[0]
        if len(clean) > best_len:
            best_len = len(clean)
            best_frame = frame
    return best_frame


def _translate_nt(nt: str) -> str:
    """Translate a NT string in-frame (frame 1 of the given string)."""
    from abverse._codons import CODON_TABLE
    aa_parts: list[str] = []
    for i in range(0, len(nt) - 2, 3):
        codon = nt[i : i + 3].upper()
        aa_parts.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa_parts)


# ── FASTA helpers ─────────────────────────────────────────────────────────────

def _read_nt_fasta(fasta_path: str) -> dict[str, str]:
    """Return {gene_name: nt_sequence} from a FASTA file."""
    seqs: dict[str, str] = {}
    current_id: Optional[str] = None
    current_parts: list[str] = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_parts)
                current_id = line[1:].split()[0]
                current_parts = []
            else:
                current_parts.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(current_parts)
    return seqs


def _sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# ── MMseqs2 database builder ──────────────────────────────────────────────────

def _make_mmseqs_protein_db(fasta_path: str, db_prefix: str) -> None:
    """Run `mmseqs createdb` to build a protein database."""
    import subprocess
    from abutils.bin import get_path as get_binary_path

    mmseqs_bin = get_binary_path("mmseqs")
    cmd = [mmseqs_bin, "createdb", fasta_path, db_prefix, "--dbtype", "1"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"mmseqs createdb failed:\n{result.stderr}"
        )


# ── AA FASTA builder ──────────────────────────────────────────────────────────

def _build_v_aa_fasta(
    nt_seqs: dict[str, str],
    output_fasta: str,
) -> dict[str, int]:
    """Translate V genes (always frame 1). Returns {gene_name: frame}.

    Sequences are trimmed at the first stop codon so that pseudogene alleles
    with internal stop codons produce a valid (truncated) AA entry rather than
    a FASTA record with embedded '*' characters.
    """
    frame_map: dict[str, int] = {}
    with open(output_fasta, "w") as fh:
        for gene, nt in nt_seqs.items():
            aa = _translate_nt(nt)
            # Trim at first stop codon (handles pseudogene alleles)
            aa = aa.split("*")[0]
            if not aa:
                continue  # entire sequence is a stop — skip
            fh.write(f">{gene}\n{aa}\n")
            frame_map[gene] = 1
    return frame_map


def _build_j_aa_fasta(
    nt_seqs: dict[str, str],
    output_fasta: str,
) -> dict[str, int]:
    """Translate J genes with automatic frame detection. Returns {gene_name: frame}."""
    frame_map: dict[str, int] = {}
    with open(output_fasta, "w") as fh:
        for gene, nt in nt_seqs.items():
            frame = _detect_j_frame(nt)
            offset = frame - 1
            aa = _translate_nt(nt[offset:])
            aa = aa.split("*")[0]  # trim at first stop codon
            fh.write(f">{gene}\n{aa}\n")
            frame_map[gene] = frame
    return frame_map


# ── Cache management ──────────────────────────────────────────────────────────

def _db_dir(species: str, receptor: str, segment: str) -> str:
    return os.path.join(_CACHE_ROOT, receptor, species, f"{segment}_aa")


def _is_cache_valid(db_dir: str, source_checksum: str) -> bool:
    checksum_path = os.path.join(db_dir, "checksum.sha256")
    if not os.path.isfile(checksum_path):
        return False
    with open(checksum_path) as fh:
        stored = fh.read().strip()
    # Check that the MMseqs2 DB files exist too
    db_prefix = os.path.join(db_dir, "mmseqs_db")
    if not os.path.isfile(db_prefix):
        return False
    return stored == source_checksum


def _write_cache_artifacts(
    db_dir: str,
    aa_fasta: str,
    frame_map: dict[str, int],
    source_checksum: str,
    db_prefix: str,
) -> None:
    os.makedirs(db_dir, exist_ok=True)
    _make_mmseqs_protein_db(aa_fasta, db_prefix)
    with open(os.path.join(db_dir, "frame_map.json"), "w") as fh:
        json.dump(frame_map, fh)
    with open(os.path.join(db_dir, "checksum.sha256"), "w") as fh:
        fh.write(source_checksum)


# ── Locate source NT germlines from abstar ────────────────────────────────────

def _abstar_ungapped_path(receptor: str, species: str) -> str:
    import abstar
    base = os.path.dirname(abstar.__file__)
    return os.path.join(base, "germline_dbs", receptor, species, "ungapped")


# ── Public API ────────────────────────────────────────────────────────────────

def get_germline_aa_db_path(species: str, receptor: str, segment: str) -> str:
    """Return path to the MMseqs2 protein DB prefix for the given segment."""
    return os.path.join(_db_dir(species, receptor, segment), "mmseqs_db")


def build_germline_aa_db(
    species: str = "human",
    receptor: str = "bcr",
    force_rebuild: bool = False,
) -> dict:
    """Build (or load from cache) AA germline databases for V and J genes.

    Returns a dict with keys:
        v_db_path, j_db_path     — MMseqs2 DB prefixes
        v_nt_seqs, j_nt_seqs    — {gene_name: nt_sequence}
        v_frame_map, j_frame_map — {gene_name: reading_frame_int}
    """
    ungapped_dir = _abstar_ungapped_path(receptor, species)
    v_fasta_path = os.path.join(ungapped_dir, "v.fasta")
    j_fasta_path = os.path.join(ungapped_dir, "j.fasta")

    if not os.path.isfile(v_fasta_path):
        raise FileNotFoundError(f"V germline FASTA not found: {v_fasta_path}")
    if not os.path.isfile(j_fasta_path):
        raise FileNotFoundError(f"J germline FASTA not found: {j_fasta_path}")

    v_nt_seqs = _read_nt_fasta(v_fasta_path)
    j_nt_seqs = _read_nt_fasta(j_fasta_path)

    result: dict = {
        "v_nt_seqs": v_nt_seqs,
        "j_nt_seqs": j_nt_seqs,
    }

    for segment, nt_seqs, source_fasta, build_fn in [
        ("v", v_nt_seqs, v_fasta_path, _build_v_aa_fasta),
        ("j", j_nt_seqs, j_fasta_path, _build_j_aa_fasta),
    ]:
        db_dir = _db_dir(species, receptor, segment)
        db_prefix = os.path.join(db_dir, "mmseqs_db")
        checksum = _sha256_file(source_fasta)

        if not force_rebuild and _is_cache_valid(db_dir, checksum):
            with open(os.path.join(db_dir, "frame_map.json")) as fh:
                frame_map = json.load(fh)
        else:
            os.makedirs(db_dir, exist_ok=True)
            aa_fasta = os.path.join(db_dir, "sequences.fasta")
            frame_map = build_fn(nt_seqs, aa_fasta)
            _write_cache_artifacts(db_dir, aa_fasta, frame_map, checksum, db_prefix)

        result[f"{segment}_db_path"] = db_prefix
        result[f"{segment}_frame_map"] = frame_map

    return result
