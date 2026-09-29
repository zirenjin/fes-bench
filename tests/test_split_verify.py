from pathlib import Path

from fes_bench.splits.verify import verify


def test_frozen_splits_have_valid_payload_hashes_and_disjoint_frames():
    root = Path(__file__).resolve().parents[1]
    report = verify(root / "data" / "processed" / "splits")
    assert report["status"] == "pass"
    assert set(report["splits"]) == {"temp_extrap", "phase_lopo", "system_loso"}
