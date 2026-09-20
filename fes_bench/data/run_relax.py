"""Compatibility launcher for representative relaxation on PyTorch 2.6+."""

from __future__ import annotations

import torch

# e3nn ships this constant file with the installed environment. PyTorch 2.6
# changed torch.load's default to weights_only=True; `slice` is the sole
# built-in required by the trusted e3nn constants archive on thu-GenSi.
torch.serialization.add_safe_globals([slice])

from .relax_representatives import main


if __name__ == "__main__":
    raise SystemExit(main())
