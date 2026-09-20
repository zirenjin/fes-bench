from __future__ import annotations

import csv
import json
from pathlib import Path

from fes_bench.qh.materialize import FQH_FIELDS, materialize


def test_materialize_preserves_diagnostic_status_and_canonical_columns(tmp_path: Path) -> None:
    data = tmp_path / "data" / "toy" / "alpha"
    source = tmp_path / "source"
    data.mkdir(parents=True)
    source.mkdir()
    (data / "structure.extxyz").write_text(
        '1\nLattice="2 0 0 0 3 0 0 0 4" Properties=species:S:1:pos:R:3\nX 0 0 0\n', encoding="utf-8"
    )
    (data / "reference_G.csv").write_text(
        "T_K,P_GPa,G_eV_per_atom,level,source\n300,1.5,-1,reference,test\n", encoding="utf-8"
    )
    (source / "alpha_fqh.csv").write_text(
        "T_K,F_QH_eV_per_atom,E_static_eV_per_atom,volume_scale,min_frequency_THz\n300,-2.5,-2.0,1.02,-0.2\n",
        encoding="utf-8",
    )
    (source / "qh_summary.json").write_text(
        json.dumps({"checkpoint": "model.pt", "phases": {"alpha": {"minimum_frequency_THz_by_volume": {"1.02": -0.2}}}}),
        encoding="utf-8",
    )
    config = tmp_path / "config.yaml"
    config.write_text(json.dumps({"data_root": "data", "source_root": "source", "system": "toy", "phases": ["alpha"], "qh_reliable": {"alpha": False}}), encoding="utf-8")

    materialize(config)

    with (data / "fqh.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
        assert tuple(rows[0]) == FQH_FIELDS
        assert float(rows[0]["P_GPa"]) == 1.5
        assert float(rows[0]["V_min_A3_per_atom"]) == 24.48
        assert float(rows[0]["F_vib_eV_per_atom"]) == -0.5
    report = json.loads((data / "phonon_report.json").read_text(encoding="utf-8"))
    assert report["qh_reliable"] is False
    assert report["diagnostic_only"] is True
    assert report["imaginary_mode_fraction"] is None
