"""Extract one-shot DPA E0/volume inputs for Bartel, never trains."""

from pathlib import Path

from ase.io import read
import torch

torch.serialization.add_safe_globals([slice])

from deepmd.calculator import DP


def main() -> None:
    root = Path("data")
    calculator = DP(
        model="external/checkpoints/DPA-3.1-3M.pt",
        head="Domains_Alloy",
    )
    for phase in ("hcp", "bcc"):
        atoms = read(root / "hf" / phase / "structure.extxyz")
        atoms.calc = calculator
        print(phase, len(atoms), atoms.get_potential_energy() / len(atoms), atoms.get_volume() / len(atoms), flush=True)


if __name__ == "__main__":
    main()
