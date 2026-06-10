"""Tests for _errors.py — the ReverseTranslationError public exception."""

import pytest

from abverse import ReverseTranslationError


class TestReverseTranslationError:
    def test_exposes_structured_failures(self):
        failures = [
            {"seq_id": "MR-72_HC", "kind": "invalid_residue",
             "detail": "'O' at position 40"},
        ]
        err = ReverseTranslationError(failures)
        assert err.failures == failures

    def test_is_an_exception(self):
        err = ReverseTranslationError([
            {"seq_id": "x", "kind": "invalid_residue", "detail": "'O' at position 1"},
        ])
        assert isinstance(err, Exception)
        with pytest.raises(ReverseTranslationError):
            raise err

    def test_renders_invalid_residue_message(self):
        failures = [
            {"seq_id": "MR-72_HC", "kind": "invalid_residue",
             "detail": "'O' at position 40"},
            {"seq_id": "MR-72_LC", "kind": "invalid_residue",
             "detail": "'O' at positions 4, 89, 105"},
        ]
        msg = str(ReverseTranslationError(failures))
        assert "2 of" in msg
        assert "invalid residues" in msg
        assert "only the 20 standard amino acids" in msg
        assert "MR-72_HC: 'O' at position 40" in msg
        assert "MR-72_LC: 'O' at positions 4, 89, 105" in msg

    def test_renders_reconstruction_error_message(self):
        failures = [
            {"seq_id": "seq_0", "kind": "reconstruction_error",
             "detail": "translate(output) != input_aa"},
        ]
        msg = str(ReverseTranslationError(failures))
        assert "failed reconstruction" in msg
        assert "seq_0: translate(output) != input_aa" in msg
