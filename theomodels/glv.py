"""
theomodels/glv.py

Generalized Lotka-Volterra (GLV) model with multiplicative noise.

SDE (Itô form):
    dX_i = X_i (r_i + sum_j A_ij X_j) dt + X_i sum_k Sigma_ik dW_k

Integration scheme: Milstein (second-order correction for multiplicative noise).

The Milstein correction for species i with multiplicative diffusion is:
    0.5 * X_i^2 * sum_k Sigma_ik^2 * dt

It exposes two public objects:

    GLVParams      dataclass holding all model parameters as tensors
    GLVSimulator   class that builds parameters from a Hydra DictConfig
                   and runs the Milstein integrator

The only interface with the rest of the project is:

    simulator = GLVSimulator(cfg)
    t, X = simulator.integrate()

where t has shape [T] and X has shape [T, n_species], both on the device
specified in cfg.model.device.

Hydra config keys consumed from cfg.model:
    n_species   int
    r           list[float]          intrinsic growth rates
    self_reg    float                diagonal of A (self-regulation, must be < 0)
    device      str                  "cpu" or "cuda"
    seed        int

Hydra config keys consumed from cfg.experiment:
    interaction_matrix   list[list[float]] | null
        Full A matrix. If null, only self-regulation is used (without_interactions).
        For mean_field, all off-diagonal entries should be equal.
        For with_interactions, entries are drawn externally and passed here.
    sigma_matrix         list[list[float]]
        Full n x n diffusion matrix Sigma. Diagonal entries control the
        per-species noise intensity. Off-diagonal entries introduce correlated
        environmental fluctuations across species.
    x0                   list[float]          initial abundances

Hydra config keys consumed from cfg.simulation:
    t_start   float
    t_end     float
    dt        float
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import torch
from omegaconf import DictConfig


# ---------------------------------------------------------------------------
# Parameter container
# ---------------------------------------------------------------------------

@dataclass
class GLVParams:
    """
    All GLV model parameters stored as PyTorch tensors on a single device.

    Attributes
    ----------
    r : Tensor, shape [n]
        Intrinsic growth rates.
    A : Tensor, shape [n, n]
        Interaction matrix. A[i, j] is the effect of species j on species i.
        The diagonal A[i, i] must be strictly negative (self-regulation).
    Sigma : Tensor, shape [n, n]
        Diffusion matrix. The instantaneous covariance of the noise per unit
        time is Sigma @ Sigma.T. A diagonal Sigma means independent
        environmental fluctuations; off-diagonal entries couple species.
    x0 : Tensor, shape [n]
        Initial abundances. All entries must be strictly positive.
    n : int
        Number of species, derived from r.
    device : torch.device
        Device on which all tensors reside.
    """

    r: torch.Tensor
    A: torch.Tensor
    Sigma: torch.Tensor
    x0: torch.Tensor

    @property
    def n(self) -> int:
        return self.r.shape[0]

    @property
    def device(self) -> torch.device:
        return self.r.device

    def validate(self) -> None:
        """
        Check tensor shapes and ecological constraints.
        Raises ValueError with an informative message on any violation.
        """
        n = self.n

        if self.A.shape != (n, n):
            raise ValueError(
                f"A must be shape ({n}, {n}), got {tuple(self.A.shape)}"
            )
        if self.Sigma.shape != (n, n):
            raise ValueError(
                f"Sigma must be shape ({n}, {n}), got {tuple(self.Sigma.shape)}"
            )
        if self.x0.shape != (n,):
            raise ValueError(
                f"x0 must be shape ({n},), got {tuple(self.x0.shape)}"
            )

        # Self-regulation: all diagonal entries must be strictly negative
        diag = torch.diag(self.A)
        if not torch.all(diag < 0):
            bad = (diag >= 0).nonzero(as_tuple=True)[0].tolist()
            raise ValueError(
                f"A[i,i] must be < 0 for all i. Violated at species: {bad}"
            )

        # Initial abundances must be positive (log-normal noise assumption)
        if not torch.all(self.x0 > 0):
            bad = (self.x0 <= 0).nonzero(as_tuple=True)[0].tolist()
            raise ValueError(
                f"x0 must be > 0 for all species. Violated at indices: {bad}"
            )


# ---------------------------------------------------------------------------
# ODE drift  (deterministic part, reused by Milstein)
# ---------------------------------------------------------------------------

def _glv_drift(X: torch.Tensor, params: GLVParams) -> torch.Tensor:
    """
    Compute the deterministic drift f(X) of the GLV model.

    f_i(X) = X_i * (r_i + sum_j A_ij X_j)

    Parameters
    ----------
    X : Tensor, shape [n]
        Current species abundances. Must be non-negative.
    params : GLVParams

    Returns
    -------
    Tensor, shape [n]
        Drift vector f(X).
    """
    # interaction term: shape [n]
    interaction = params.A @ X            # [n, n] x [n] -> [n]
    growth_rate = params.r + interaction  # [n]
    return X * growth_rate               # element-wise: [n]


# ---------------------------------------------------------------------------
# Milstein integrator
# ---------------------------------------------------------------------------

def _milstein_step(
    X: torch.Tensor,
    params: GLVParams,
    dt: float,
    sqrt_dt: float,
    generator: torch.Generator,
) -> torch.Tensor:
    """
    Advance the state by one Milstein step.

    SDE (Itô):
        dX_i = f_i(X) dt + X_i [Sigma dW]_i

    where [Sigma dW]_i = sum_k Sigma_ik dW_k and dW ~ N(0, dt I).

    Milstein scheme:
        X_i(t+dt) = X_i(t)
                  + f_i(X) dt                             (Euler drift)
                  + X_i [Sigma dW]_i                      (diffusion)
                  + 0.5 X_i^2 sum_k Sigma_ik^2 dt        (Milstein correction)

    The correction term 0.5 X_i^2 sum_k Sigma_ik^2 dt arises from the
    Itô-Taylor expansion of the diffusion coefficient g_i(X) = X_i with
    respect to X_i. For a full derivation see Kloeden & Platen (1992),
    Chapter 10.

    After the step, abundances are clamped to a small positive value to
    prevent numerical extinction from carrying the trajectory to X_i < 0,
    which has no ecological meaning under the log-normal noise model.

    Parameters
    ----------
    X : Tensor, shape [n]
        Current state.
    params : GLVParams
    dt : float
        Time step.
    sqrt_dt : float
        Precomputed sqrt(dt) for efficiency.
    generator : torch.Generator
        Seeded random number generator for reproducibility.

    Returns
    -------
    Tensor, shape [n]
        Updated state X(t + dt).
    """
    n = params.n

    # Independent Wiener increments: dW ~ N(0, dt) for each of the n
    # driving noise dimensions
    dW = torch.zeros(n, dtype=X.dtype, device=X.device)
    dW.normal_(mean=0.0, std=sqrt_dt, generator=generator)  # shape [n]

    # Diffusion term: [Sigma dW]_i = sum_k Sigma_ik dW_k
    # Shape: [n, n] x [n] -> [n]
    sigma_dW = params.Sigma @ dW                             # [n]

    # Drift
    drift = _glv_drift(X, params)                           # [n]

    # Milstein correction: 0.5 * X_i^2 * sum_k Sigma_ik^2 * dt
    # sum_k Sigma_ik^2 is the row-wise sum of squares of Sigma
    sigma_sq_rowsum = (params.Sigma ** 2).sum(dim=1)        # [n]
    milstein_correction = 0.5 * (X ** 2) * sigma_sq_rowsum * dt  # [n]

    # Full Milstein update
    X_new = X + drift * dt + X * sigma_dW + milstein_correction

    # Clamp to avoid numerical extinction (minimum abundance = 1e-10)
    X_new = torch.clamp(X_new, min=1e-10)

    return X_new


# ---------------------------------------------------------------------------
# Public simulator class
# ---------------------------------------------------------------------------

class GLVSimulator:
    """
    Builds GLV parameters from a Hydra DictConfig and runs the Milstein
    integrator.

    Usage
    -----
    Called from simulate.py after Hydra has resolved the config:

        simulator = GLVSimulator(cfg)
        t, X = simulator.integrate()

    Parameters
    ----------
    cfg : DictConfig
        Fully resolved Hydra config. The keys consumed are documented in the
        module docstring.
    """

    def __init__(self, cfg: DictConfig) -> None:
        self.cfg = cfg
        self.params = self._build_params(cfg)
        self.params.validate()

        # Time grid
        self.t_start: float = float(cfg.simulation.t_start)
        self.t_end: float   = float(cfg.simulation.t_end)
        self.dt: float      = float(cfg.simulation.dt)
        self.n_steps: int   = math.ceil(
            (self.t_end - self.t_start) / self.dt
        )

        # Reproducible random generator — separate from the global state
        # so that other parts of the codebase are not affected
        self._generator = torch.Generator(device=self.params.device)
        self._generator.manual_seed(int(cfg.model.seed))

    # ------------------------------------------------------------------
    # Parameter construction
    # ------------------------------------------------------------------

    @staticmethod
    def _build_params(cfg: DictConfig) -> GLVParams:
        """
        Construct GLVParams from the resolved Hydra config.

        The interaction matrix A is built as follows:
          - Diagonal is always set to cfg.model.self_reg (must be < 0).
          - If cfg.experiment.interaction_matrix is provided, its
            off-diagonal entries override the diagonal-only default.
          - If null, only self-regulation is used (without_interactions).

        The diffusion matrix Sigma is taken directly from
        cfg.experiment.sigma_matrix.
        """
        device = torch.device(cfg.model.device)
        dtype  = torch.float64   # ecological simulations need float64

        n = int(cfg.model.n_species)

        # Intrinsic growth rates
        r = torch.tensor(list(cfg.model.r), dtype=dtype, device=device)

        # Interaction matrix A
        if cfg.experiment.interaction_matrix is None:
            # without_interactions: diagonal only
            A = torch.zeros(n, n, dtype=dtype, device=device)
        else:
            A = torch.tensor(
                [list(row) for row in cfg.experiment.interaction_matrix],
                dtype=dtype,
                device=device,
            )

        # Enforce self-regulation on the diagonal regardless of what was
        # passed — this prevents accidental positive diagonal values
        A.fill_diagonal_(float(cfg.model.self_reg))

        # Diffusion matrix Sigma
        Sigma = torch.tensor(
            [list(row) for row in cfg.experiment.sigma_matrix],
            dtype=dtype,
            device=device,
        )

        # Initial conditions
        x0 = torch.tensor(
            list(cfg.experiment.x0),
            dtype=dtype,
            device=device,
        )

        return GLVParams(r=r, A=A, Sigma=Sigma, x0=x0)

    # ------------------------------------------------------------------
    # Integration
    # ------------------------------------------------------------------

    def integrate(self) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Run the Milstein integrator from t_start to t_end.

        Returns
        -------
        t : Tensor, shape [T]
            Time points. T = n_steps + 1 (includes t=0).
        X : Tensor, shape [T, n_species]
            Species abundances at each time point.

        Notes
        -----
        The trajectory is stored entirely in memory. For very long
        simulations (large n_steps) or many species, consider streaming
        directly to HDF5 from simulate.py rather than returning the full
        tensor.
        """
        params   = self.params
        dt       = self.dt
        sqrt_dt  = math.sqrt(dt)
        n_steps  = self.n_steps
        device   = params.device
        dtype    = params.r.dtype

        # Pre-allocate output tensors
        T = n_steps + 1
        t_out = torch.zeros(T, dtype=dtype, device=device)
        X_out = torch.zeros(T, params.n, dtype=dtype, device=device)

        # Set initial conditions
        X_out[0] = params.x0
        t_out[0] = self.t_start

        # Milstein loop
        X = params.x0.clone()
        for step in range(n_steps):
            X = _milstein_step(X, params, dt, sqrt_dt, self._generator)
            t_out[step + 1] = self.t_start + (step + 1) * dt
            X_out[step + 1] = X

        return t_out, X_out

    # ------------------------------------------------------------------
    # Convenience: multiple realisations (ensemble)
    # ------------------------------------------------------------------

    def integrate_ensemble(
        self, n_realisations: int
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Run n_realisations independent trajectories from the same x0.

        Each realisation uses a different seed derived from cfg.model.seed
        so the ensemble is reproducible but internally diverse.

        Returns
        -------
        t : Tensor, shape [T]
            Shared time grid.
        X_ensemble : Tensor, shape [n_realisations, T, n_species]
            One trajectory per realisation.

        This is useful for computing statistics over noise realisations
        (e.g. mean trajectory, variance bands) to be stored in HDF5 under
        simulations/{model}/{experiment}/ensemble/.
        """
        params  = self.params
        dt      = self.dt
        sqrt_dt = math.sqrt(dt)
        n_steps = self.n_steps
        device  = params.device
        dtype   = params.r.dtype
        base_seed = int(self.cfg.model.seed)

        T = n_steps + 1
        X_ensemble = torch.zeros(
            n_realisations, T, params.n, dtype=dtype, device=device
        )

        # Use a shared time grid from the first realisation
        t_out = torch.zeros(T, dtype=dtype, device=device)
        t_out[0] = self.t_start
        for step in range(n_steps):
            t_out[step + 1] = self.t_start + (step + 1) * dt

        for k in range(n_realisations):
            gen = torch.Generator(device=device)
            gen.manual_seed(base_seed + k)    # distinct but reproducible seed

            X = params.x0.clone()
            X_ensemble[k, 0] = X
            for step in range(n_steps):
                X = _milstein_step(X, params, dt, sqrt_dt, gen)
                X_ensemble[k, step + 1] = X

        return t_out, X_ensemble
