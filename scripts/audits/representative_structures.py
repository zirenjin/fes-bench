"""Audit representative-structure sources and relaxations — plan experiment E8."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from common import write_audit


def relative_source(value: object) -> object:
    if isinstance(value, str) and value.startswith("/"):
        return "external/" + Path(value).name
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--repo-root", type=Path, default=Path(".")); args = parser.parse_args(); root = args.repo_root.resolve()
    rows, inputs, missing = [], [], []
    for path in sorted((root / "data").glob("*/system.json")):
        system = path.parent.name; inputs.append(path)
        for phase in json.loads(path.read_text(encoding="utf-8")).get("phases", []):
            meta_path = path.parent / phase / "meta.json"
            if not meta_path.exists(): missing.append({"field": f"{system}:{phase}", "reason": "phase meta.json missing"}); continue
            meta = json.loads(meta_path.read_text(encoding="utf-8")); inputs.append(meta_path); rep = meta.get("representative", {}); energy = rep.get("energy_eV_per_atom", {}); rows.append({"system": system, "phase": phase, "status": rep.get("status", ""), "source": relative_source(rep.get("source", rep.get("source_structure", ""))), "calculator": rep.get("calculator", ""), "head": rep.get("head", ""), "energy_before_eV_per_atom": energy.get("before", ""), "energy_after_eV_per_atom": energy.get("after", ""), "max_displacement_A": rep.get("max_displacement_A", ""), "converged": rep.get("converged", "")})
    write_audit(root, "representative_structures", ["system", "phase", "status", "source", "calculator", "head", "energy_before_eV_per_atom", "energy_after_eV_per_atom", "max_displacement_A", "converged"], rows, inputs, missing)
    return 0


if __name__ == "__main__": raise SystemExit(main())
