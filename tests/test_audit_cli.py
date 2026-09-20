from __future__ import annotations

from pathlib import Path

from fes_bench.data.audit import main


def test_empty_phase0_audit_is_a_valid_scaffold(tmp_path: Path, capsys) -> None:
    config = tmp_path / "audit.yaml"
    config.write_text('{"data_root": "' + str(tmp_path) + '", "systems": []}', encoding="utf-8")

    assert main(["--config", str(config)]) == 0
    output = capsys.readouterr().out
    assert "systems=0" in output
