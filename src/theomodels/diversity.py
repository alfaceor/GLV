"""Diversity and divergence metrics for GLV trajectory tensors.

All functions expect trajectories of shape ``(n_trials, time, n_species)``
and operate on the *species* axis (``dim=-1``).  Ensemble statistics
(mean, std, quantiles) are then reduced over the *trials* axis (``dim=0``).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import torch


# ---------------------------------------------------------------------------
# Low-level helpers (single distribution or pair)
# ---------------------------------------------------------------------------

def relative_abundance(X: torch.Tensor, eps: float = 1e-12) -> torch.Tensor:
    """Relative abundances (probability simplex) along the species axis.

    Args:
        X: Abundances with species on the last dimension. Shape (..., n_species).
        eps: Small constant to avoid division by zero.

    Returns:
        Tensor of same shape as ``X``, summing to 1 along the last dim.
    """
    if torch.any(X < 0):
        raise ValueError("Abundances must be non-negative.")
    total = X.sum(dim=-1, keepdim=True)
    return X / (total + eps)


def shannon_diversity(
    X: torch.Tensor,
    base: float | None = None,
    eps: float = 1e-12,
) -> torch.Tensor:
    """Shannon diversity H = -∑ p_i log(p_i) over relative abundances.

    Args:
        X: Abundances with species on the last dimension. Shape (..., n_species).
        base: Logarithm base.  ``None`` (default) → natural log (nats);
              ``2`` → bits.
        eps: Small constant for numerical stability.

    Returns:
        Shannon diversity, shape ``(...)`` (species dim reduced).
    """
    p = relative_abundance(X, eps=eps)
    log_p = torch.log(p + eps)
    if base is not None:
        log_p = log_p / math.log(base)
    return -(p * log_p).sum(dim=-1)


def kl_divergence(
    P: torch.Tensor,
    Q: torch.Tensor,
    base: float | None = None,
    eps: float = 1e-12,
    from_abundances: bool = True,
) -> torch.Tensor:
    """Kullback–Leibler divergence D_KL(P ‖ Q).

    Args:
        P: Distribution P — abundances or probabilities, species on last dim.
        Q: Distribution Q — must be broadcastable with ``P``.
        base: Logarithm base (``None`` → nats).
        eps: Small constant for numerical stability.
        from_abundances: Convert P and Q to relative abundances first when
            ``True``; assume they are already valid probability vectors when
            ``False``.

    Returns:
        KL divergence, shape ``(...)`` (species dim reduced).
    """
    if from_abundances:
        p = relative_abundance(P, eps=eps)
        q = relative_abundance(Q, eps=eps)
    else:
        p, q = P, Q
    log_ratio = torch.log((p + eps) / (q + eps))
    if base is not None:
        log_ratio = log_ratio / math.log(base)
    return (p * log_ratio).sum(dim=-1)


def js_divergence(
    P: torch.Tensor,
    Q: torch.Tensor,
    base: float | None = None,
    eps: float = 1e-12,
    from_abundances: bool = True,
) -> torch.Tensor:
    """Jensen–Shannon divergence JSD(P ‖ Q).

    Symmetric and bounded in ``[0, log(2)]`` (nats).

    Args:
        P: Distribution P — abundances or probabilities, species on last dim.
        Q: Distribution Q — must be broadcastable with ``P``.
        base: Logarithm base (``None`` → nats).
        eps: Small constant for numerical stability.
        from_abundances: See :func:`kl_divergence`.

    Returns:
        JS divergence, shape ``(...)`` (species dim reduced).
    """
    if from_abundances:
        p = relative_abundance(P, eps=eps)
        q = relative_abundance(Q, eps=eps)
    else:
        p, q = P, Q
    m = 0.5 * (p + q)
    return 0.5 * (
        kl_divergence(p, m, base=base, eps=eps, from_abundances=False)
        + kl_divergence(q, m, base=base, eps=eps, from_abundances=False)
    )


# ---------------------------------------------------------------------------
# Trajectory-level helpers — operate on (n_trials, time, n_species)
# ---------------------------------------------------------------------------

def traj_relative_abundance(
    traj: torch.Tensor,
    eps: float = 1e-12,
) -> torch.Tensor:
    """Relative abundance at every time-step for every realization.

    Args:
        traj: Shape ``(n_trials, time, n_species)``.
        eps: Numerical stability constant.

    Returns:
        Shape ``(n_trials, time, n_species)``.
    """
    return relative_abundance(traj, eps=eps)


def traj_shannon(
    traj: torch.Tensor,
    base: float | None = None,
    eps: float = 1e-12,
) -> torch.Tensor:
    """Shannon diversity at every time-step for every realization.

    Args:
        traj: Shape ``(n_trials, time, n_species)``.
        base: Logarithm base (``None`` → nats).
        eps: Numerical stability constant.

    Returns:
        Shape ``(n_trials, time)``.
    """
    return shannon_diversity(traj, base=base, eps=eps)


def traj_js_vs_ensemble_mean(
    traj: torch.Tensor,
    base: float | None = None,
    eps: float = 1e-12,
) -> torch.Tensor:
    """JSD between each realization and the ensemble-mean distribution.

    At every time-step the reference distribution is the mean relative
    abundance across all trials, making JSD a measure of how far each
    realization strays from the ensemble consensus.

    Args:
        traj: Shape ``(n_trials, time, n_species)``.
        base: Logarithm base (``None`` → nats).
        eps: Numerical stability constant.

    Returns:
        Shape ``(n_trials, time)``.
    """
    p = relative_abundance(traj, eps=eps)          # (n_trials, time, n_species)
    ref = p.mean(dim=0, keepdim=True)              # (1, time, n_species)
    return js_divergence(p, ref, base=base, eps=eps, from_abundances=False)


# ---------------------------------------------------------------------------
# Ensemble statistics — mean / std / quantiles across trials
# ---------------------------------------------------------------------------

@dataclass
class EnsembleStats:
    """Container for ensemble statistics of a per-realization metric.

    All tensors share shape ``(time,)`` (or ``(time, n_species)`` for
    relative-abundance statistics) and live on the same device as the
    input trajectory.

    Attributes:
        mean:     Ensemble mean across trials at each time-step.
        std:      Ensemble standard deviation across trials.
        q25:      25th percentile across trials.
        median:   50th percentile across trials.
        q75:      75th percentile across trials.
        q05:      5th percentile across trials.
        q95:      95th percentile across trials.
    """
    mean:   torch.Tensor
    std:    torch.Tensor
    q05:    torch.Tensor
    q25:    torch.Tensor
    median: torch.Tensor
    q75:    torch.Tensor
    q95:    torch.Tensor

    def to_dict(self) -> dict[str, torch.Tensor]:
        return {
            "mean":   self.mean,
            "std":    self.std,
            "q05":    self.q05,
            "q25":    self.q25,
            "median": self.median,
            "q75":    self.q75,
            "q95":    self.q95,
        }


def _ensemble_stats(metric: torch.Tensor) -> EnsembleStats:
    """Compute ensemble statistics over the first (trials) dimension.

    Args:
        metric: Tensor of shape ``(n_trials, ...)``.

    Returns:
        :class:`EnsembleStats` with each field having shape ``(...)``.
    """
    qs = torch.quantile(
        metric.float(),
        torch.tensor([0.05, 0.25, 0.50, 0.75, 0.95], device=metric.device),
        dim=0,
    )
    return EnsembleStats(
        mean=metric.mean(dim=0),
        std=metric.std(dim=0),
        q05=qs[0],
        q25=qs[1],
        median=qs[2],
        q75=qs[3],
        q95=qs[4],
    )


# ---------------------------------------------------------------------------
# Main entry point — compute everything in one call
# ---------------------------------------------------------------------------

@dataclass
class DiversityResults:
    """All per-realization diversity metrics and their ensemble statistics.

    Per-realization tensors
    -----------------------
    rel_abundance : ``(n_trials, time, n_species)``
    shannon       : ``(n_trials, time)``
    jsd           : ``(n_trials, time)``

    Ensemble statistics (over the trials axis)
    ------------------------------------------
    rel_abundance_stats : :class:`EnsembleStats` with shape
                          ``(time, n_species)`` per field.
    shannon_stats       : :class:`EnsembleStats` with shape ``(time,)``
                          per field.
    jsd_stats           : :class:`EnsembleStats` with shape ``(time,)``
                          per field.
    """
    # Per-realization
    rel_abundance: torch.Tensor
    shannon:       torch.Tensor
    jsd:           torch.Tensor
    # Ensemble
    rel_abundance_stats: EnsembleStats
    shannon_stats:       EnsembleStats
    jsd_stats:           EnsembleStats


def compute_diversity(
    traj: torch.Tensor,
    base: float | None = None,
    eps: float = 1e-12,
) -> DiversityResults:
    """Compute diversity metrics and ensemble statistics from a trajectory tensor.

    Args:
        traj: Simulation output of shape ``(n_trials, time, n_species)``.
        base: Logarithm base for Shannon and JSD (``None`` → nats, ``2`` →
              bits).
        eps:  Numerical stability constant forwarded to all metric functions.

    Returns:
        :class:`DiversityResults` containing per-realization tensors and their
        ensemble statistics (mean, std, 5th/25th/50th/75th/95th quantiles).

    Example::

        traj = simulate(A, r, X0, dt=0.01, steps=1000, n_trials=50)
        results = compute_diversity(traj)
        print(results.shannon_stats.mean.shape)   # (1001,)
        print(results.rel_abundance.shape)        # (50, 1001, n_species)
    """
    rel_ab = traj_relative_abundance(traj, eps=eps)
    shannon = traj_shannon(traj, base=base, eps=eps)
    jsd = traj_js_vs_ensemble_mean(traj, base=base, eps=eps)

    return DiversityResults(
        rel_abundance=rel_ab,
        shannon=shannon,
        jsd=jsd,
        rel_abundance_stats=_ensemble_stats(rel_ab),
        shannon_stats=_ensemble_stats(shannon),
        jsd_stats=_ensemble_stats(jsd),
    )