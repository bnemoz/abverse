# Public exception for reverse-translation failures.

from __future__ import annotations

__all__ = ["ReverseTranslationError"]


class ReverseTranslationError(Exception):
    """Raised when one or more sequences cannot be reverse-translated.

    Attributes
    ----------
    failures : list[dict]
        One entry per failed sequence, each with keys:
          - seq_id : str
          - kind   : "invalid_residue" | "reconstruction_error"
          - detail : str   (human-readable, e.g. "'O' at positions 4, 89, 105")
    """

    def __init__(self, failures: list[dict]):
        self.failures = failures
        super().__init__(self._render())

    def _render(self) -> str:
        invalid = [f for f in self.failures if f["kind"] == "invalid_residue"]
        recon = [f for f in self.failures if f["kind"] == "reconstruction_error"]
        other = [
            f for f in self.failures
            if f["kind"] not in ("invalid_residue", "reconstruction_error")
        ]

        blocks: list[str] = []
        total = len(self.failures)

        if invalid:
            header = (
                f"{len(invalid)} of {total} input sequences contain invalid "
                f"residues\n(only the 20 standard amino acids are accepted):"
            )
            lines = [f"  {f['seq_id']}: {f['detail']}" for f in invalid]
            blocks.append("\n".join([header, *lines]))

        if recon:
            header = f"{len(recon)} sequence(s) failed reconstruction:"
            lines = [f"  {f['seq_id']}: {f['detail']}" for f in recon]
            blocks.append("\n".join([header, *lines]))

        if other:
            lines = [f"  {f['seq_id']}: {f['detail']}" for f in other]
            blocks.append("\n".join(["Reverse translation failed:", *lines]))

        return "\n\n".join(blocks) if blocks else "Reverse translation failed."
