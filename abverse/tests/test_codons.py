"""Tests for _codons.py — codon lookup table and selection logic."""

import pytest

from abverse._codons import (
    CODON_TABLE,
    HUMAN_PREFERRED_CODON,
    OPTIMAL_CODON_TABLE,
    _SYNONYMOUS,
    _hamming,
    fallback_codon,
    optimal_codon,
)


class TestCodonTable:
    def test_all_64_codons_present(self):
        assert len(CODON_TABLE) == 64

    def test_known_translations(self):
        assert CODON_TABLE["ATG"] == "M"
        assert CODON_TABLE["TGG"] == "W"
        assert CODON_TABLE["TAA"] == "*"
        assert CODON_TABLE["TAG"] == "*"
        assert CODON_TABLE["TGA"] == "*"

    def test_synonymous_coverage(self):
        # Every non-stop AA should appear in _SYNONYMOUS
        aas = set(CODON_TABLE.values()) - {"*"}
        assert aas == set(_SYNONYMOUS.keys()) - {"*"}

    def test_synonymous_round_trips(self):
        for aa, codons in _SYNONYMOUS.items():
            for codon in codons:
                assert CODON_TABLE[codon] == aa


class TestHamming:
    def test_identical(self):
        assert _hamming("ATG", "ATG") == 0

    def test_one_diff(self):
        assert _hamming("ATG", "TTG") == 1

    def test_all_diff(self):
        assert _hamming("ATG", "CAC") == 3


class TestOptimalCodonTable:
    def test_size(self):
        # 20 AAs × 64 germline codons = 1280 entries
        assert len(OPTIMAL_CODON_TABLE) == 1280

    def test_output_encodes_target_aa(self):
        for (target_aa, germ_codon), result_codon in OPTIMAL_CODON_TABLE.items():
            assert CODON_TABLE.get(result_codon) == target_aa, (
                f"optimal_codon({target_aa!r}, {germ_codon!r}) = {result_codon!r} "
                f"does not encode {target_aa!r}"
            )

    def test_hamming_minimised(self):
        """For each entry, no synonymous codon should have a lower Hamming distance."""
        from abverse._codons import _HUMAN_FREQ

        for (target_aa, germ_codon), result_codon in OPTIMAL_CODON_TABLE.items():
            best_hamming = _hamming(result_codon, germ_codon)
            for synonym in _SYNONYMOUS[target_aa]:
                assert _hamming(synonym, germ_codon) >= best_hamming, (
                    f"Found shorter Hamming for ({target_aa}, {germ_codon}): "
                    f"{synonym} ({_hamming(synonym, germ_codon)}) < {result_codon} ({best_hamming})"
                )

    def test_tie_broken_by_frequency(self):
        # CTG and CTC both encode L and differ from CTG germline by 0 and 1 — CTG wins
        result = OPTIMAL_CODON_TABLE[("L", "CTG")]
        assert result == "CTG"


class TestOptimalCodonFunction:
    def test_zero_hamming_wins(self):
        # If germline codon already encodes target AA, it should be returned as-is
        # ATG encodes M; optimal_codon('M', 'ATG') == 'ATG'
        assert optimal_codon("M", "ATG") == "ATG"

    def test_nonstandard_aa_returns_nnn(self):
        assert optimal_codon("X", "ATG") == "NNN"
        assert optimal_codon("B", "ATG") == "NNN"
        assert optimal_codon("Z", "ATG") == "NNN"

    def test_ambiguous_germline_codon(self):
        # Should not crash; falls back to human preferred codon
        result = optimal_codon("A", "NNN")
        assert CODON_TABLE.get(result) == "A"

    def test_wrong_length_germline_codon(self):
        result = optimal_codon("A", "GC")
        assert CODON_TABLE.get(result) == "A"


class TestFallbackCodon:
    def test_preferred_codons(self):
        for aa, codon in HUMAN_PREFERRED_CODON.items():
            if aa == "*":
                continue
            assert fallback_codon(aa) == codon

    def test_nonstandard_aa(self):
        assert fallback_codon("X") == "NNN"
        assert fallback_codon("B") == "NNN"
