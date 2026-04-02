"""Generalized Lotka-Volterra simulation using PyTorch and Hydra."""

from __future__ import annotations

from dataclasses import dataclass

import hydra
from hydra.core.config_store import ConfigStore

import torch
import matplotlib.pyplot as plt


def simulate(A: torch.Tensor, r: torch.Tensor, X0: torch.Tensor, dt: float, steps: int) -> torch.Tensor:
    """Euler-step simulation of dX_i = X_i (r_i + sum_j A_ij X_j)."""
    X = X0.clone()
    traj = torch.zeros(steps + 1, X.size(0), dtype=X.dtype, device=X.device)
    traj[0] = X
    for i in range(steps):
        dX = X * (r + A @ X)
        X = X + dt * dX
        traj[i + 1] = X
    return traj


def random_system(n_species: int, device: torch.device | None = None, dtype: torch.dtype = torch.float32):
    """Generate random interaction matrix, growth rates, and initial state."""
    device = device or torch.device("cpu")
    A = torch.randn(n_species, n_species, device=device, dtype=dtype)
    r = torch.randn(n_species, device=device, dtype=dtype)
    X0 = torch.rand(n_species, device=device, dtype=dtype)
    return A, r, X0


def plot_trajectories(traj: torch.Tensor) -> None:
    """Plot the trajectories of all species."""
    tt = torch.arange(traj.size(0))  # Time points
    plt.plot(tt.cpu().numpy(), traj.cpu().numpy())
    plt.xlabel("Time")
    plt.ylabel("Abundance")
    plt.title("Generalized Lotka-Volterra Simulation")
    plt.show()

@dataclass
class Config:
    n_species: int = 3
    dt: float = 0.01
    steps: int = 100


cs = ConfigStore.instance()
cs.store(name="glv_config", node=Config)

@hydra.main(version_base=None, config_name="glv_config")
def main(cfg: Config) -> None:
    print("Configuration:", cfg)
    A, r, X0 = random_system(cfg.n_species)
    # tt = torch.arange(cfg.steps + 1) * cfg.dt
    traj = simulate(A, r, X0, cfg.dt, cfg.steps)
    plot_trajectories(traj)
    print("Final state:", traj[-1])


if __name__ == "__main__":
    main()

