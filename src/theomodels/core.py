"""Core helper utilities for GLV simulations."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from datetime import datetime

import torch
import matplotlib.pyplot as plt
import h5py

import numpy as np

def simulate(
    A: torch.Tensor,
    r: torch.Tensor,
    X0: torch.Tensor,
    dt: float,
    steps: int,
    sigma: torch.Tensor | float | None = None,
    n_trials: int = 1,
    f_t: torch.Tensor | None = None,
) -> torch.Tensor:
    """Euler-Maruyama simulation of dX_i = X_i (r_i + sum_j A_ij X_j) with multiplicative noise."""
    if steps < 0:
        raise ValueError("steps must be non-negative")
    if dt <= 0:
        raise ValueError("dt must be positive")
    if f_t is not None:
        if f_t.shape[0] != steps or f_t.shape[1] != r.shape[0]:
            raise ValueError(
                f"Invalid external for f_t tensor shape: {f_t.shape}",
                f"Expected tensor shape {(steps, r.shape[0])}"
            )

    X = X0.clone()
    n_species = X.size(0)
    traj = torch.zeros(n_trials, steps + 1, n_species, dtype=X.dtype, device=X.device)
    # Broadcast of initial conditions
    traj[:, 0] = X

    r = r.to(device=X.device, dtype=X.dtype)
    A = A.to(device=X.device, dtype=X.dtype)

    if sigma is None:
        noise_strength = None
    elif isinstance(sigma, torch.Tensor):
        noise_strength = sigma.to(device=X.device, dtype=X.dtype)
    else:
        noise_strength = torch.full((n_species,), float(sigma), device=X.device, dtype=X.dtype)

    sqrt_dt = math.sqrt(dt)
    for trial in range(n_trials):
        for i in range(steps):
            drift = X * (r + A @ X)
            if f_t is not None:
                drift += f_t[i]

            if noise_strength is None:
                X = X + dt * drift
            else:
                stochastic = noise_strength * X * sqrt_dt * torch.randn_like(X)
                X = X + dt * drift + stochastic
            X = torch.clamp(X, min=1e-6)
            traj[trial, i + 1] = X
    return traj

def b_is_feasible_Q(A, r) -> bool:
    """Is the system with interaction matrix A and growth rate r feasible, meaning $A^{-1} r  > 0$

    Args:
        A (_type_): interaction matrix
        r (_type_): growth rate

    Returns:
        bool: _description_
    """
    A_inv = np.linalg.inv(A)
    n_star = np.matmul(A_inv, r)
    return all(-n_star > 0)
    

def random_gaussian_system(
    n_species: int,
    device: torch.device | None = None,
    dtype: torch.dtype = torch.float32,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Sample Gaussian A and r without guarantees on feasibility or stability."""
    device = device or torch.device("cpu")
    A = torch.randn(n_species, n_species, device=device, dtype=dtype)
    r = torch.randn(n_species, device=device, dtype=dtype)
    X0 = torch.rand(n_species, device=device, dtype=dtype)
    return A, r, X0

# TODO: create a function that generates a random A and r. Based on a defined x_star that will be a stationary state
def random_feasible_system_x_star(
    x_star: torch.Tensor,
    device: torch.device | None = None,
    dtype: torch.dtype = torch.float32
) -> tuple[torch.Tensor, torch.Tensor]:
    """Generate a random species interaction and growth rate based on the stationary state x_star 

    Args:
        x_star (torch.Tensor): Unidimensional tensor final stationary state for system.
        device (torch.device | None, optional): _description_. Defaults to None.
        dtype (torch.dtype, optional): _description_. Defaults to torch.float32.

    Returns:
        tuple[torch.Tensor, torch.Tensor]: Interaction matrix, growth rate for each species.
    """
    n_species = x_star.shape[0]
    device = device or torch.device("cpu")

    B = torch.randn(n_species, n_species, device=device, dtype=dtype)
    A = -(B @ B.T) - 0.1 * torch.eye(n_species, device=device, dtype=dtype)
    r = -(A @ x_star)
    return A, r



def random_feasible_system(
    n_species: int,
    device: torch.device | None = None,
    dtype: torch.dtype = torch.float32,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Generate A and r with a positive, linearly stable equilibrium."""
    device = device or torch.device("cpu")

    x_star = 0.5 + torch.rand(n_species, device=device, dtype=dtype)

    B = torch.randn(n_species, n_species, device=device, dtype=dtype)
    A = -(B @ B.T) - 0.1 * torch.eye(n_species, device=device, dtype=dtype)

    r = -(A @ x_star)

    X0 = x_star * (0.8 + 0.4 * torch.rand(n_species, device=device, dtype=dtype))

    return A, r, X0

def random_system(
    n_species: int,
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

    return generator(n_species=n_species, device=device, dtype=dtype)

