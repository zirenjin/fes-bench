"""Two-layer temperature plus phase-ID MLP comparison baseline."""

from __future__ import annotations

import numpy as np
import torch


class PhaseIdMLP:
    """Fit a scalar free-energy curve from only temperature and phase identity.

    This deliberately has no structure, composition, pressure, or QH input.
    It is a temperature-function comparison baseline, not a transferable
    physical free-energy model.
    """

    def __init__(
        self,
        phases: list[str],
        seed: int = 11,
        steps: int = 1000,
        learning_rate: float = 1e-3,
    ) -> None:
        if not phases or len(set(phases)) != len(phases):
            raise ValueError("phases must be a non-empty list of unique names")
        if steps <= 0 or learning_rate <= 0.0:
            raise ValueError("steps and learning_rate must be positive")
        self.phases = list(phases)
        self.steps = int(steps)
        self.learning_rate = float(learning_rate)
        # Construction must not perturb a caller's global random stream.
        with torch.random.fork_rng(devices=[]):
            torch.manual_seed(seed)
            self.net = torch.nn.Sequential(
                torch.nn.Linear(1 + len(self.phases), 32),
                torch.nn.Tanh(),
                torch.nn.Linear(32, 32),
                torch.nn.Tanh(),
                torch.nn.Linear(32, 1),
            )

    def _x(self, temperature_K: np.ndarray, phases: list[str]) -> torch.Tensor:
        temperature = np.asarray(temperature_K, dtype=float).reshape(-1)
        if len(temperature) != len(phases) or not np.isfinite(temperature).all():
            raise ValueError("temperature and phases must have equal finite lengths")
        try:
            phase_index = [self.phases.index(phase) for phase in phases]
        except ValueError as exc:
            raise ValueError("encountered phase absent from the configured phase vocabulary") from exc
        one_hot = np.zeros((len(phases), len(self.phases)), dtype=np.float32)
        one_hot[np.arange(len(phases)), phase_index] = 1.0
        return torch.tensor(np.column_stack((temperature / 1000.0, one_hot)), dtype=torch.float32)

    def fit(self, temperature_K: np.ndarray, phases: list[str], energy_eV_per_atom: np.ndarray) -> None:
        """Fit on frozen training points only."""

        x = self._x(temperature_K, phases)
        energy = np.asarray(energy_eV_per_atom, dtype=float).reshape(-1)
        if len(energy) != len(phases) or not np.isfinite(energy).all():
            raise ValueError("energy and phases must have equal finite lengths")
        y = torch.tensor(energy, dtype=torch.float32).reshape(-1, 1)
        optimizer = torch.optim.Adam(self.net.parameters(), lr=self.learning_rate)
        self.net.train()
        for _ in range(self.steps):
            optimizer.zero_grad()
            loss = torch.nn.functional.mse_loss(self.net(x), y)
            loss.backward()
            optimizer.step()

    def predict(self, temperature_K: np.ndarray, phases: list[str]) -> np.ndarray:
        """Return eV/atom predictions for the supplied phase-labelled points."""

        self.net.eval()
        with torch.no_grad():
            return self.net(self._x(temperature_K, phases)).squeeze(1).numpy()
