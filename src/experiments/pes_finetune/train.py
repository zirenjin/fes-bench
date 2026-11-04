"""Launch a DFT-only DPA PES fine-tune with hard free-energy guards — plan PES step 3."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def assert_free_energy_blind(data_root: Path) -> None:
    forbidden = []
    for path in data_root.rglob("*"):
        if not path.is_file():
            continue
        name = path.name.lower()
        if name == "reference_g.csv" or "thermodynamic" in name or name in {"fqh.csv", "qh_summary.json"}:
            forbidden.append(str(path))
    if forbidden:
        raise AssertionError("PES fine-tuning input contains forbidden free-energy/thermodynamic files: " + ", ".join(forbidden[:10]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--head", required=True, choices=("Domains_SSE_PBE", "Domains_Alloy"))
    parser.add_argument("--system", required=True, choices=("sio2", "hf", "ti", "zr"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dp", default="dp")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    data_root = args.data_root.resolve()
    config = args.config.resolve()
    checkpoint = args.checkpoint.resolve()
    if not data_root.exists():
        raise FileNotFoundError(data_root)
    assert_free_energy_blind(data_root)
    config_text = config.read_text(encoding="utf-8")
    if any(token in config_text.lower() for token in ("reference_g.csv", "fqh.csv", "thermodynamic")):
        raise AssertionError("training config contains a forbidden free-energy input")
    if not checkpoint.is_file():
        raise FileNotFoundError(checkpoint)
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    provenance = {
        "status": "dry_run" if args.dry_run else "launching",
        "free_energy_blind": True,
        "system": args.system,
        "functional": "PBE-D3(BJ)" if args.system == "sio2" else "PBE",
        "head": args.head,
        "checkpoint_path": str(checkpoint),
        "checkpoint_sha256": sha256(checkpoint),
        "data_root": str(data_root),
        "config_path": str(config),
        "config_sha256": sha256(config),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "forbidden_inputs_checked": ["reference_G.csv", "fqh.csv", "thermodynamic tables"],
        "dataloader_num_workers": 0,
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    if args.dry_run:
        print(output / "provenance.json")
        return 0
    command = [args.dp, "train", str(config)]
    with (output / "train.log").open("w", encoding="utf-8") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    provenance["status"] = "complete" if completed.returncode == 0 else "failed"
    provenance["returncode"] = completed.returncode
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
