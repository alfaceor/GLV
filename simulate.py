"""Generalized Lotka-Volterra simulation using PyTorch and Hydra."""

from __future__ import annotations

from dataclasses import dataclass
import math
from datetime import datetime
from pathlib import Path

import hydra
from hydra.core.config_store import ConfigStore

import torch
import matplotlib.pyplot as plt
import h5py


def simulate(
    A: torch.Tensor,
    r: torch.Tensor,
    X0: torch.Tensor,
    dt: float,
    steps: int,
    sigma: torch.Tensor | float | None = None,
) -> torch.Tensor:
    """Euler-Maruyama simulation of dX_i = X_i (r_i + sum_j A_ij X_j) with multiplicative noise."""
    if steps < 0:
        raise ValueError("steps must be non-negative")
    if dt <= 0:
        raise ValueError("dt must be positive")

    X = X0.clone()
    n_species = X.size(0)
    traj = torch.zeros(steps + 1, n_species, dtype=X.dtype, device=X.device)
    traj[0] = X

    r = r.to(device=X.device, dtype=X.dtype)
    A = A.to(device=X.device, dtype=X.dtype)

    if sigma is None:
        noise_strength = None
    elif isinstance(sigma, torch.Tensor):
        noise_strength = sigma.to(device=X.device, dtype=X.dtype)
    else:
        noise_strength = torch.full((n_species,), float(sigma), device=X.device, dtype=X.dtype)

    sqrt_dt = math.sqrt(dt)

    for i in range(steps):
        drift = X * (r + A @ X)
        if noise_strength is None:
            X = X + dt * drift
        else:
            stochastic = noise_strength * X * sqrt_dt * torch.randn_like(X)
            X = X + dt * drift + stochastic
        traj[i + 1] = X
    return traj


def random_gaussian_system(
    n_species: int,
    noise_std: float = 0.1,
    device: torch.device | None = None,
    dtype: torch.dtype = torch.float32,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Sample Gaussian A and r without guarantees on feasibility or stability."""
    device = device or torch.device("cpu")
    A = torch.randn(n_species, n_species, device=device, dtype=dtype)
    r = torch.randn(n_species, device=device, dtype=dtype)
    X0 = torch.rand(n_species, device=device, dtype=dtype)
    sigma = torch.full((n_species,), float(noise_std), device=device, dtype=dtype)
    return A, r, X0, sigma


def random_feasible_system(
    n_species: int,
    noise_std: float = 0.1,
    device: torch.device | None = None,
    dtype: torch.dtype = torch.float32,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Generate A and r with a positive, linearly stable equilibrium."""
    device = device or torch.device("cpu")

    # Pick a strictly positive equilibrium
    x_star = 0.5 + torch.rand(n_species, device=device, dtype=dtype)

    # Build a symmetric negative definite interaction matrix
    B = torch.randn(n_species, n_species, device=device, dtype=dtype)
    A = -(B @ B.T) - 0.1 * torch.eye(n_species, device=device, dtype=dtype)

    # Set growth rates so that x_star is an equilibrium: r = -A x_star
    r = -(A @ x_star)

    # Sample initial condition near the equilibrium
    X0 = x_star * (0.8 + 0.4 * torch.rand(n_species, device=device, dtype=dtype))

    sigma = torch.full((n_species,), float(noise_std), device=device, dtype=dtype)
    return A, r, X0, sigma


def random_system(
    n_species: int,
    noise_std: float = 0.1,
    device: torch.device | None = None,
    dtype: torch.dtype = torch.float32,
    mode: str = "feasible",
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Dispatch between Gaussian and feasible/stable random systems."""
    generators = {
        "gaussian": random_gaussian_system,
        "feasible": random_feasible_system,
    }
    try:
        generator = generators[mode]
    except KeyError as exc:  # pragma: no cover - simple input guard
        raise ValueError(f"unknown random system mode '{mode}'") from exc

    return generator(n_species=n_species, noise_std=noise_std, device=device, dtype=dtype)



def plot_trajectories(
    traj: torch.Tensor,
    ax: plt.Axes,
    linestyle: str = "-",
    label: str | None = None,
) -> None:
    """Plot the trajectories of all species."""
    tt = torch.arange(traj.size(0))  # Time points
    ax.plot(tt.cpu().numpy(), traj.cpu().numpy(), linestyle=linestyle, label=label)
    ax.set_xlabel("Time")
    ax.set_ylabel("Abundance")
    ax.set_title("Generalized Lotka-Volterra Simulation")


def plot_matrix_A(A: torch.Tensor, ax: plt.Axes) -> None:
    """Plot the interaction matrix."""
    im = ax.imshow(A.cpu().numpy(), cmap="coolwarm")
    ax.set_title("Interaction matrix A")
    ax.set_xlabel("Species j")
    ax.set_ylabel("Species i")
    plt.colorbar(im, ax=ax, orientation="vertical", shrink=0.8)


@dataclass
class Config:
    n_species: int = 64
    dt: float = 0.01
    steps: int = 50
    noise_std: float = 0.1
    seed: int | None = 42
    system_mode: str = "feasible"
    device: str = "auto"
    save_results: bool = True
    output_dir: str = "outputs"
    run_name: str = "glv_sim"




def save_to_hdf5(
    filepath: Path,
    cfg: Config,
    deterministic: torch.Tensor,
    stochastic: torch.Tensor,
    A: torch.Tensor,
    r: torch.Tensor,
    X0: torch.Tensor,
    sigma: torch.Tensor,
    metadata: dict[str, str],
) -> None:
    """Persist simulation outputs and configuration to an HDF5 file."""
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(filepath, "w") as h5:
        h5.create_dataset("trajectories/deterministic", data=deterministic.cpu().numpy())
        h5.create_dataset("trajectories/stochastic", data=stochastic.cpu().numpy())
        h5.create_dataset("system/A", data=A.cpu().numpy())
        h5.create_dataset("system/r", data=r.cpu().numpy())
        h5.create_dataset("system/X0", data=X0.cpu().numpy())
        h5.create_dataset("system/sigma", data=sigma.cpu().numpy())

        cfg_group = h5.create_group("config")
        for key, value in vars(cfg).items():
            if value is None:
                continue
            cfg_group.attrs[key] = str(value)

        meta_group = h5.create_group("metadata")
        for key, value in metadata.items():
            meta_group.attrs[key] = value


cs = ConfigStore.instance()
cs.store(name="glv_config", node=Config)
@hydra.main(version_base=None, config_name="glv_config")
# @hydra.main(version_base=None, config_path="simulconfigs", config_name="config")
def main(cfg: Config) -> None:
    if cfg.seed is not None:
        torch.manual_seed(cfg.seed)

    if cfg.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(cfg.device)

    print("Configuration:", cfg)
    print(f"Using device: {device}")

    A, r, X0, sigma = random_system(
        cfg.n_species,
        noise_std=cfg.noise_std,
        device=device,
        mode=cfg.system_mode,
    )
    # TODO: For each simulation:
    # 1. Run the deterministic baseline first.
    # 2. Run n_trials stochastic trials; save the average trajectory.
    # 3. Save the trials with maximum and minimum deviation from the deterministic baseline.
    deter_traj = simulate(A, r, X0, cfg.dt, cfg.steps)
    stoch_traj = simulate(A, r, X0, cfg.dt, cfg.steps, sigma=sigma)

    if cfg.save_results:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path(cfg.output_dir)
        filename = output_dir / f"{cfg.run_name}_{timestamp}.h5"
        metadata = {
            "device": str(device),
            "timestamp": timestamp,
            "system_mode": cfg.system_mode,
        }
        save_to_hdf5(
            filename,
            cfg,
            deterministic=deter_traj,
            stochastic=stoch_traj,
            A=A,
            r=r,
            X0=X0,
            sigma=sigma,
            metadata=metadata,
        )
        print(f"Saved simulation results to {filename}")

    fig, ax = plt.subplots(2, 1, figsize=(10, 8))
    plot_trajectories(deter_traj, ax=ax[0], linestyle="--", label="deterministic")
    plot_trajectories(stoch_traj, ax=ax[0], linestyle="-", label="stochastic")
    plot_matrix_A(A, ax=ax[1])

    print("Final state (deterministic):", deter_traj[-1].cpu())
    print("Final state (stochastic):", stoch_traj[-1].cpu())

    ax[0].legend()
    plt.show()


if __name__ == "__main__":
    main()

