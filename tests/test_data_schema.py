from __future__ import annotations

import json
from pathlib import Path

import pytest

from fes_bench.data import load
from fes_bench.data.schema import SchemaError


def _write_valid_phase(root: Path) -> None:
    phase_dir = root / "toy" / "alpha"
    phase_dir.mkdir(parents=True)
    (root / "toy" / "system.json").write_text(
        json.dumps(
            {
                "phases": ["alpha"],
                "type_map": ["X"],
                "reference_grid": {"T_K": [300.0], "P_GPa": [0.0]},
                "truth_level": "test",
            }
        ),
        encoding="utf-8",
    )
    (phase_dir / "structure.extxyz").write_text("1\nProperties=species:S:1:pos:R:3\nX 0 0 0\n", encoding="utf-8")
    (phase_dir / "reference_G.csv").write_text(
        "T_K,P_GPa,G_eV_per_atom,level,source\n300,0,-1.25,test,unit-test\n",
        encoding="utf-8",
    )
    (phase_dir / "meta.json").write_text(
        json.dumps(
            {
                "functional": "test",
                "dispersion": "none",
                "supercell": "1x1x1",
                "kpoints": "gamma",
                "convergence": "test",
                "size_error_eV_per_atom": 0.0,
                "method": "TI",
                "doi": "test",
                "notes": "fixture",
            }
        ),
        encoding="utf-8",
    )


def test_load_returns_normalized_phase(tmp_path: Path) -> None:
    _write_valid_phase(tmp_path)

    phase = load("toy", "alpha", tmp_path)

    assert phase.structure.name == "structure.extxyz"
    assert phase.G_table[0].T_K == 300.0
    assert phase.G_table[0].G_eV_per_atom == -1.25
    assert phase.meta["method"] == "TI"


def test_load_rejects_noncanonical_reference_header(tmp_path: Path) -> None:
    _write_valid_phase(tmp_path)
    (tmp_path / "toy" / "alpha" / "reference_G.csv").write_text(
        "T_K,G_eV_per_atom\n300,-1.25\n", encoding="utf-8"
    )

    with pytest.raises(SchemaError, match="exactly these columns"):
        load("toy", "alpha", tmp_path)
