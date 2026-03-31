"""
tests/theomodels/test_glv.py

Unit tests for theomodels/glv.py.

Run with:
    pytest tests/theomodels/test_glv.py -v

Dependencies:
    pytest
    torch
    omegaconf

No Hydra runtime is needed — tests build DictConfig objects directly using
OmegaConf.create(), which is what Hydra resolves to internally. This keeps
tests fast and free of filesystem side effects.
"""

import math

import pytest
import torch

from omegaconf import OmegaConf

from theomodels.glv import GLVParams, GLVSimulator, _glv_drift, _milstein_step


# ---------------------------------------------------------------------------
# Shared config factories
# ---------------------------------------------------------------------------
# These functions build minimal DictConfig objects that match exactly the
# keys GLVSimulator reads from cfg. They are kept as plain functions rather
# than pytest fixtures so they can be called with arguments inside tests.


def _base_cfg(
    n: int = 2,
    r: list[float] | None = None,
    self_reg: float = -0.5,
    seed: int = 0,
    device: str = "cpu",
    x0: list[float] | None = None,
    interaction_matrix=None,
    sigma_matrix: list[list[float]] | None = None,
    t_start: float = 0.0,
    t_end: float = 5.0,
    dt: float = 0.01,
):
    """
    Build a minimal OmegaConf DictConfig that GLVSimulator can consume.
    All arguments have sensible defaults for a 2-species system.
    """
    if r is None:
        r = [1.0] * n
    if x0 is None:
        x0 = [1.0] * n
    if sigma_matrix is None:
        # Small diagonal noise by default
        sigma_matrix = [[0.1 if i == j else 0.0 for j in range(n)]
                        for i in range(n)]

    return OmegaConf.create({
        "model": {
            "n_species": n,
            "r": r,
            "self_reg": self_reg,
            "seed": seed,
            "device": device,
        },
        "experiment": {
            "interaction_matrix": interaction_matrix,
            "sigma_matrix": sigma_matrix,
            "x0": x0,
        },
        "simulation": {
            "t_start": t_start,
            "t_end": t_end,
            "dt": dt,
        },
    })


# ---------------------------------------------------------------------------
# GLVParams — construction and validation
# ---------------------------------------------------------------------------

class TestGLVParams:

    def test_n_property_matches_r_length(self):
        r     = torch.tensor([1.0, 0.5, -0.2])
        A     = -0.1 * torch.eye(3)
        Sigma = 0.05 * torch.eye(3)
        x0    = torch.tensor([1.0, 2.0, 0.5])
        params = GLVParams(r=r, A=A, Sigma=Sigma, x0=x0)
        assert params.n == 3

    def test_device_property_returns_tensor_device(self):
        r     = torch.tensor([1.0])
        A     = torch.tensor([[-0.5]])
        Sigma = torch.tensor([[0.1]])
        x0    = torch.tensor([1.0])
        params = GLVParams(r=r, A=A, Sigma=Sigma, x0=x0)
        assert params.device == torch.device("cpu")

    def test_validate_passes_for_valid_params(self):
        n = 3
        r     = torch.ones(n)
        A     = -0.5 * torch.eye(n)
        Sigma = 0.1 * torch.eye(n)
        x0    = torch.ones(n)
        params = GLVParams(r=r, A=A, Sigma=Sigma, x0=x0)
        params.validate()   # should not raise

    def test_validate_raises_for_wrong_A_shape(self):
        r     = torch.ones(3)
        A     = torch.eye(2) * -0.5          # wrong shape
        Sigma = torch.eye(3) * 0.1
        x0    = torch.ones(3)
        params = GLVParams(r=r, A=A, Sigma=Sigma, x0=x0)
        with pytest.raises(ValueError, match="A must be shape"):
            params.validate()

    def test_validate_raises_for_wrong_Sigma_shape(self):
        n = 3
        r     = torch.ones(n)
        A     = -0.5 * torch.eye(n)
        Sigma = torch.eye(2) * 0.1           # wrong shape
        x0    = torch.ones(n)
        params = GLVParams(r=r, A=A, Sigma=Sigma, x0=x0)
        with pytest.raises(ValueError, match="Sigma must be shape"):
            params.validate()

    def test_validate_raises_for_non_negative_diagonal_A(self):
        n = 3
        r     = torch.ones(n)
        A     = torch.zeros(n, n)
        A[1, 1] = 0.1                        # positive diagonal — invalid
        Sigma = 0.1 * torch.eye(n)
        x0    = torch.ones(n)
        params = GLVParams(r=r, A=A, Sigma=Sigma, x0=x0)
        with pytest.raises(ValueError, match="A\[i,i\] must be < 0"):
            params.validate()

    def test_validate_raises_for_non_positive_x0(self):
        n = 2
        r     = torch.ones(n)
        A     = -0.5 * torch.eye(n)
        Sigma = 0.1 * torch.eye(n)
        x0    = torch.tensor([1.0, -0.5])    # negative abundance — invalid
        params = GLVParams(r=r, A=A, Sigma=Sigma, x0=x0)
        with pytest.raises(ValueError, match="x0 must be > 0"):
            params.validate()


# ---------------------------------------------------------------------------
# _glv_drift — deterministic vector field
# ---------------------------------------------------------------------------

class TestGLVDrift:

    def _make_params(self, n, r_vals, A_vals, device="cpu"):
        dtype = torch.float64
        r     = torch.tensor(r_vals, dtype=dtype, device=device)
        A     = torch.tensor(A_vals, dtype=dtype, device=device)
        Sigma = 0.1 * torch.eye(n, dtype=dtype, device=device)
        x0    = torch.ones(n, dtype=dtype, device=device)
        return GLVParams(r=r, A=A, Sigma=Sigma, x0=x0)

    def test_without_interactions_equals_logistic(self):
        """
        With A diagonal only, each species grows logistically.
        f_i = X_i * (r_i + A_ii * X_i).
        """
        n   = 3
        r   = [1.0, 0.5, 0.8]
        A   = [[-0.5, 0.0, 0.0],
               [0.0, -0.3, 0.0],
               [0.0,  0.0, -0.4]]
        params = self._make_params(n, r, A)
        X = torch.tensor([1.0, 2.0, 0.5], dtype=torch.float64)

        drift = _glv_drift(X, params)

        expected = torch.tensor([
            X[0] * (r[0] + A[0][0] * X[0]),
            X[1] * (r[1] + A[1][1] * X[1]),
            X[2] * (r[2] + A[2][2] * X[2]),
        ], dtype=torch.float64)

        assert torch.allclose(drift, expected, atol=1e-12)

    def test_drift_zero_at_carrying_capacity_single_species(self):
        """
        For a single species without interactions:
        carrying capacity K = -r / A_11.
        At X = K, f(X) = X * (r + A*X) = X * (r - r) = 0.
        """
        r_val  = 1.0
        A_val  = -0.5
        K      = -r_val / A_val   # = 2.0

        params = self._make_params(
            n=1,
            r_vals=[r_val],
            A_vals=[[A_val]],
        )
        X = torch.tensor([K], dtype=torch.float64)
        drift = _glv_drift(X, params)

        assert torch.allclose(drift, torch.zeros(1, dtype=torch.float64), atol=1e-12)

    def test_drift_is_zero_at_X_zero(self):
        """
        GLV drift is always zero when species abundance is zero
        (multiplicative structure: f_i = X_i * (...) = 0 when X_i = 0).
        """
        n   = 2
        r   = [1.0, -0.5]
        A   = [[-0.5, 0.1], [-0.2, -0.3]]
        params = self._make_params(n, r, A)
        X = torch.zeros(n, dtype=torch.float64)

        drift = _glv_drift(X, params)
        assert torch.allclose(drift, torch.zeros(n, dtype=torch.float64), atol=1e-12)

    def test_drift_shape(self):
        n = 5
        r = [1.0] * n
        A = [[-0.5 if i == j else 0.05 for j in range(n)] for i in range(n)]
        params = self._make_params(n, r, A)
        X = torch.ones(n, dtype=torch.float64)
        assert _glv_drift(X, params).shape == (n,)


# ---------------------------------------------------------------------------
# _milstein_step — single integration step
# ---------------------------------------------------------------------------

class TestMilsteinStep:

    def _make_params_and_gen(self, n=2, sigma_scale=0.0, seed=42):
        """Helper: build params with controllable noise level."""
        dtype = torch.float64
        r     = torch.ones(n, dtype=dtype)
        A     = -0.5 * torch.eye(n, dtype=dtype)
        Sigma = sigma_scale * torch.eye(n, dtype=dtype)
        x0    = torch.ones(n, dtype=dtype)
        params = GLVParams(r=r, A=A, Sigma=Sigma, x0=x0)
        gen = torch.Generator()
        gen.manual_seed(seed)
        return params, gen

    def test_zero_noise_equals_euler_forward(self):
        """
        With Sigma = 0, the Milstein step reduces to an Euler step:
        X_new = X + f(X) * dt.
        Milstein correction vanishes (0.5 * X^2 * 0 * dt = 0).
        Diffusion term vanishes (X * 0 * dW = 0).
        """
        n = 2
        params, gen = self._make_params_and_gen(n=n, sigma_scale=0.0)
        X  = torch.tensor([1.5, 0.8], dtype=torch.float64)
        dt = 0.01

        X_new = _milstein_step(X, params, dt, math.sqrt(dt), gen)
        X_euler = X + _glv_drift(X, params) * dt

        assert torch.allclose(X_new, X_euler, atol=1e-12)

    def test_output_shape_unchanged(self):
        n = 4
        params, gen = self._make_params_and_gen(n=n, sigma_scale=0.05)
        X  = torch.ones(n, dtype=torch.float64)
        dt = 0.01
        X_new = _milstein_step(X, params, dt, math.sqrt(dt), gen)
        assert X_new.shape == (n,)

    def test_output_is_positive(self):
        """
        After clamping, all abundances must remain positive regardless of
        the noise realisation.
        """
        n = 3
        params, gen = self._make_params_and_gen(n=n, sigma_scale=5.0, seed=99)
        X  = torch.tensor([0.01, 0.01, 0.01], dtype=torch.float64)
        dt = 0.1
        for _ in range(50):
            X = _milstein_step(X, params, dt, math.sqrt(dt), gen)
        assert torch.all(X > 0)

    def test_different_seeds_give_different_results(self):
        """Two generators with different seeds must produce different steps."""
        n = 2
        params, gen1 = self._make_params_and_gen(n=n, sigma_scale=0.1, seed=0)
        _,       gen2 = self._make_params_and_gen(n=n, sigma_scale=0.1, seed=1)
        X  = torch.ones(n, dtype=torch.float64)
        dt = 0.01
        X1 = _milstein_step(X.clone(), params, dt, math.sqrt(dt), gen1)
        X2 = _milstein_step(X.clone(), params, dt, math.sqrt(dt), gen2)
        assert not torch.allclose(X1, X2)

    def test_same_seed_gives_identical_results(self):
        """Reproducibility: same seed must reproduce the exact same step."""
        n = 2
        params, gen1 = self._make_params_and_gen(n=n, sigma_scale=0.1, seed=7)
        _,       gen2 = self._make_params_and_gen(n=n, sigma_scale=0.1, seed=7)
        X  = torch.ones(n, dtype=torch.float64)
        dt = 0.01
        X1 = _milstein_step(X.clone(), params, dt, math.sqrt(dt), gen1)
        X2 = _milstein_step(X.clone(), params, dt, math.sqrt(dt), gen2)
        assert torch.allclose(X1, X2)


# ---------------------------------------------------------------------------
# GLVSimulator — construction from config
# ---------------------------------------------------------------------------

class TestGLVSimulatorConstruction:

    def test_builds_without_interactions(self):
        cfg = _base_cfg(n=2, interaction_matrix=None)
        sim = GLVSimulator(cfg)
        # Diagonal must be self_reg, off-diagonal must be zero
        A = sim.params.A
        assert torch.allclose(torch.diag(A),
                               torch.full((2,), -0.5, dtype=A.dtype))
        off_diag = A - torch.diag(torch.diag(A))
        assert torch.allclose(off_diag, torch.zeros_like(off_diag))

    def test_builds_with_interaction_matrix(self):
        A_cfg = [[-0.5, 0.1], [-0.2, -0.5]]
        cfg   = _base_cfg(n=2, interaction_matrix=A_cfg)
        sim   = GLVSimulator(cfg)
        # Off-diagonal entries must match config
        assert math.isclose(sim.params.A[0, 1].item(),  0.1, rel_tol=1e-9)
        assert math.isclose(sim.params.A[1, 0].item(), -0.2, rel_tol=1e-9)

    def test_self_reg_always_overrides_diagonal(self):
        """
        Even if the user passes a diagonal value in interaction_matrix,
        GLVSimulator must overwrite it with cfg.model.self_reg.
        """
        A_cfg = [[-9999.0, 0.1], [-0.2, -9999.0]]   # wrong diagonal
        cfg   = _base_cfg(n=2, interaction_matrix=A_cfg, self_reg=-0.5)
        sim   = GLVSimulator(cfg)
        assert torch.allclose(
            torch.diag(sim.params.A),
            torch.full((2,), -0.5, dtype=sim.params.A.dtype),
        )

    def test_sigma_matrix_loaded_correctly(self):
        sigma = [[0.2, 0.05], [0.05, 0.3]]
        cfg   = _base_cfg(n=2, sigma_matrix=sigma)
        sim   = GLVSimulator(cfg)
        expected = torch.tensor(sigma, dtype=torch.float64)
        assert torch.allclose(sim.params.Sigma, expected)

    def test_n_steps_computed_correctly(self):
        cfg = _base_cfg(t_start=0.0, t_end=10.0, dt=0.1)
        sim = GLVSimulator(cfg)
        assert sim.n_steps == math.ceil(10.0 / 0.1)

    def test_validate_called_on_bad_config(self):
        """Positive self_reg must raise during construction via validate()."""
        cfg = _base_cfg(n=2, self_reg=0.1)   # positive self_reg is invalid
        with pytest.raises(ValueError, match="A\[i,i\] must be < 0"):
            GLVSimulator(cfg)


# ---------------------------------------------------------------------------
# GLVSimulator.integrate — trajectory properties
# ---------------------------------------------------------------------------

class TestGLVSimulatorIntegrate:

    def test_output_shapes(self):
        n   = 3
        cfg = _base_cfg(n=n, t_start=0.0, t_end=2.0, dt=0.01,
                        x0=[1.0, 0.5, 2.0])
        sim = GLVSimulator(cfg)
        t, X = sim.integrate()
        expected_T = sim.n_steps + 1
        assert t.shape == (expected_T,)
        assert X.shape == (expected_T, n)

    def test_time_axis_starts_and_ends_correctly(self):
        cfg = _base_cfg(t_start=1.0, t_end=4.0, dt=0.05)
        sim = GLVSimulator(cfg)
        t, _ = sim.integrate()
        assert math.isclose(t[0].item(), 1.0, rel_tol=1e-9)
        # Last time point should be within dt of t_end
        assert t[-1].item() <= 4.0 + 1e-9

    def test_initial_condition_preserved(self):
        x0  = [1.5, 0.8]
        cfg = _base_cfg(n=2, x0=x0)
        sim = GLVSimulator(cfg)
        _, X = sim.integrate()
        expected = torch.tensor(x0, dtype=torch.float64)
        assert torch.allclose(X[0], expected, atol=1e-12)

    def test_all_abundances_positive(self):
        """
        Clamping must ensure X > 0 throughout the trajectory.
        Tested with large noise to stress the positivity constraint.
        """
        n     = 3
        sigma = [[0.5 if i == j else 0.0 for j in range(n)]
                 for i in range(n)]
        cfg   = _base_cfg(n=n, sigma_matrix=sigma, t_end=10.0, dt=0.01,
                          seed=123)
        sim   = GLVSimulator(cfg)
        _, X  = sim.integrate()
        assert torch.all(X > 0)

    def test_zero_noise_trajectory_is_reproducible(self):
        """Without noise the trajectory is deterministic — same every run."""
        sigma = [[0.0, 0.0], [0.0, 0.0]]
        cfg   = _base_cfg(n=2, sigma_matrix=sigma)
        t1, X1 = GLVSimulator(cfg).integrate()
        t2, X2 = GLVSimulator(cfg).integrate()
        assert torch.allclose(X1, X2)

    def test_same_seed_reproduces_stochastic_trajectory(self):
        cfg    = _base_cfg(n=2, seed=42)
        t1, X1 = GLVSimulator(cfg).integrate()
        t2, X2 = GLVSimulator(cfg).integrate()
        assert torch.allclose(X1, X2)

    def test_different_seed_gives_different_trajectory(self):
        cfg1 = _base_cfg(n=2, seed=0)
        cfg2 = _base_cfg(n=2, seed=1)
        _, X1 = GLVSimulator(cfg1).integrate()
        _, X2 = GLVSimulator(cfg2).integrate()
        assert not torch.allclose(X1, X2)

    def test_without_interactions_converges_to_carrying_capacity(self):
        """
        With zero noise and no inter-species interactions, each species
        converges to K_i = -r_i / A_ii as t -> inf.
        """
        r_vals   = [1.0, 0.5]
        self_reg = -0.5
        K_exact  = [-r / self_reg for r in r_vals]   # [2.0, 1.0]

        sigma = [[0.0, 0.0], [0.0, 0.0]]
        cfg   = _base_cfg(
            n=2,
            r=r_vals,
            self_reg=self_reg,
            sigma_matrix=sigma,
            interaction_matrix=None,
            x0=[0.1, 0.1],
            t_end=30.0,
            dt=0.01,
        )
        sim  = GLVSimulator(cfg)
        _, X = sim.integrate()

        K_tensor = torch.tensor(K_exact, dtype=torch.float64)
        assert torch.allclose(X[-1], K_tensor, atol=1e-2)

    def test_mean_field_trajectory_shape_and_positivity(self):
        n     = 4
        a_off = -0.05   # weak competition
        A_cfg = [[(-0.5 if i == j else a_off) for j in range(n)]
                 for i in range(n)]
        cfg   = _base_cfg(n=n, interaction_matrix=A_cfg, t_end=5.0)
        sim   = GLVSimulator(cfg)
        t, X  = sim.integrate()
        assert X.shape[1] == n
        assert torch.all(X > 0)


# ---------------------------------------------------------------------------
# GLVSimulator.integrate_ensemble — ensemble properties
# ---------------------------------------------------------------------------

class TestGLVSimulatorEnsemble:

    def test_ensemble_output_shapes(self):
        n   = 2
        cfg = _base_cfg(n=n, t_end=2.0, dt=0.05)
        sim = GLVSimulator(cfg)
        k   = 5
        t, X_ens = sim.integrate_ensemble(k)
        expected_T = sim.n_steps + 1
        assert t.shape == (expected_T,)
        assert X_ens.shape == (k, expected_T, n)

    def test_ensemble_realisations_differ(self):
        """Each realisation uses a different seed so they must diverge."""
        cfg      = _base_cfg(n=2, t_end=5.0)
        sim      = GLVSimulator(cfg)
        _, X_ens = sim.integrate_ensemble(n_realisations=3)
        # Realisations 0 and 1 must not be identical
        assert not torch.allclose(X_ens[0], X_ens[1])

    def test_ensemble_all_positive(self):
        n     = 3
        sigma = [[0.3 if i == j else 0.0 for j in range(n)]
                 for i in range(n)]
        cfg      = _base_cfg(n=n, sigma_matrix=sigma, t_end=5.0, seed=7)
        sim      = GLVSimulator(cfg)
        _, X_ens = sim.integrate_ensemble(n_realisations=10)
        assert torch.all(X_ens > 0)

    def test_ensemble_shared_initial_condition(self):
        """All realisations must start from the same x0."""
        x0  = [2.0, 0.5]
        cfg = _base_cfg(n=2, x0=x0)
        sim = GLVSimulator(cfg)
        _, X_ens = sim.integrate_ensemble(n_realisations=4)
        x0_tensor = torch.tensor(x0, dtype=torch.float64)
        for k in range(4):
            assert torch.allclose(X_ens[k, 0], x0_tensor, atol=1e-12)

    def test_ensemble_reproducible_across_calls(self):
        """
        Two calls to integrate_ensemble with the same cfg must give
        identical results because seeds are derived deterministically.
        """
        cfg = _base_cfg(n=2, seed=99)
        sim = GLVSimulator(cfg)
        _, X1 = sim.integrate_ensemble(n_realisations=3)
        _, X2 = sim.integrate_ensemble(n_realisations=3)
        assert torch.allclose(X1, X2)
