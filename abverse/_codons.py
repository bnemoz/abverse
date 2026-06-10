# Codon lookup table, optimal germline-informed codon selection, and fallback codon utilities.
#
# The OPTIMAL_CODON_TABLE is pre-computed at import time: for every (target_aa, germline_codon)
# pair (1280 entries = 20 AAs × 64 germline codons), it stores the synonymous codon that
# minimises Hamming distance to the germline codon, with ties broken by human codon frequency.

from __future__ import annotations

__all__ = ["optimal_codon", "fallback_codon"]

# ── Standard genetic code ────────────────────────────────────────────────────

CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*", "TGA": "*",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C",
    "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Inverted: amino acid → list of synonymous codons
_SYNONYMOUS: dict[str, list[str]] = {}
for _codon, _aa in CODON_TABLE.items():
    _SYNONYMOUS.setdefault(_aa, []).append(_codon)

# ── Human codon frequency (Homo sapiens, high-expression table) ──────────────
# Source: Kazusa codon usage database (human, RefSeq mRNA), most-frequent codon per AA.
# DNAchisel h_sapiens table used to derive these; hardcoded here to avoid runtime dependency.
HUMAN_PREFERRED_CODON: dict[str, str] = {
    "A": "GCC",
    "C": "TGC",
    "D": "GAC",
    "E": "GAG",
    "F": "TTC",
    "G": "GGC",
    "H": "CAC",
    "I": "ATC",
    "K": "AAG",
    "L": "CTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCC",
    "Q": "CAG",
    "R": "AGG",
    "S": "AGC",
    "T": "ACC",
    "V": "GTG",
    "W": "TGG",
    "Y": "TAC",
    "*": "TGA",
}

# Human codon frequency table (relative frequency, approximate).
# Used for tie-breaking in optimal_codon and ordering candidates.
_HUMAN_FREQ: dict[str, float] = {
    # Phe
    "TTT": 0.45, "TTC": 0.55,
    # Leu
    "TTA": 0.07, "TTG": 0.13, "CTT": 0.13, "CTC": 0.20, "CTA": 0.07, "CTG": 0.40,
    # Ile
    "ATT": 0.36, "ATC": 0.48, "ATA": 0.16,
    # Met
    "ATG": 1.00,
    # Val
    "GTT": 0.18, "GTC": 0.24, "GTA": 0.11, "GTG": 0.47,
    # Ser
    "TCT": 0.15, "TCC": 0.22, "TCA": 0.15, "TCG": 0.06, "AGT": 0.15, "AGC": 0.27,
    # Pro
    "CCT": 0.28, "CCC": 0.33, "CCA": 0.27, "CCG": 0.11,
    # Thr
    "ACT": 0.25, "ACC": 0.36, "ACA": 0.28, "ACG": 0.11,
    # Ala
    "GCT": 0.26, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11,
    # Tyr
    "TAT": 0.43, "TAC": 0.57,
    # Stops
    "TAA": 0.28, "TAG": 0.20, "TGA": 0.52,
    # His
    "CAT": 0.41, "CAC": 0.59,
    # Gln
    "CAA": 0.25, "CAG": 0.75,
    # Asn
    "AAT": 0.46, "AAC": 0.54,
    # Lys
    "AAA": 0.42, "AAG": 0.58,
    # Asp
    "GAT": 0.46, "GAC": 0.54,
    # Glu
    "GAA": 0.42, "GAG": 0.58,
    # Cys
    "TGT": 0.45, "TGC": 0.55,
    # Trp
    "TGG": 1.00,
    # Arg
    "CGT": 0.08, "CGC": 0.19, "CGA": 0.11, "CGG": 0.21, "AGA": 0.20, "AGG": 0.21,
    # Gly
    "GGT": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25,
}


# ── Hamming distance between two equal-length strings ────────────────────────

def _hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


# ── Pre-computed optimal codon lookup (1280 entries) ─────────────────────────
# Key: (target_aa, germline_codon)  Value: best synonymous codon

def _build_optimal_codon_lookup() -> dict[tuple[str, str], str]:
    table: dict[tuple[str, str], str] = {}
    for aa, synonyms in _SYNONYMOUS.items():
        if aa == "*":
            continue
        for germline_codon in CODON_TABLE:
            # Sort candidates by (hamming_to_germline ASC, human_freq DESC)
            best = min(
                synonyms,
                key=lambda c: (_hamming(c, germline_codon), -_HUMAN_FREQ.get(c, 0.0)),
            )
            table[(aa, germline_codon)] = best
    return table


OPTIMAL_CODON_TABLE: dict[tuple[str, str], str] = _build_optimal_codon_lookup()


# ── Public API ────────────────────────────────────────────────────────────────

def optimal_codon(target_aa: str, germline_codon: str) -> str:
    """Return the synonymous codon for *target_aa* closest (Hamming) to *germline_codon*.

    Ties broken by human codon frequency. O(1) lookup. Input is assumed
    pre-validated to the 20 standard amino acids; a non-standard residue
    raises (via :func:`fallback_codon`).
    """
    if len(germline_codon) != 3:
        return fallback_codon(target_aa)
    germline_codon = germline_codon.upper()
    result = OPTIMAL_CODON_TABLE.get((target_aa, germline_codon))
    if result is None:
        # germline_codon contains ambiguous bases — fall back to frequency
        return fallback_codon(target_aa)
    return result


def fallback_codon(target_aa: str) -> str:
    """Return the most frequent human codon for *target_aa*.

    Used for CDR3, overhangs, and positions without a valid germline codon.
    Input is assumed pre-validated; a non-standard residue raises ``ValueError``
    rather than silently producing junk.
    """
    try:
        return HUMAN_PREFERRED_CODON[target_aa]
    except KeyError:
        raise ValueError(f"No codon for non-standard residue {target_aa!r}")
