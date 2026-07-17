import pytest
import numpy as np
import torch
from theomodels.core import (
    random_feasible_system_x_star,
    b_is_feasible_Q,
    simulate
)

from theomodels.core import (  # adjust import to your actual module name
    relative_abundance,
    shannon_diversity,
    kl_divergence,
    js_divergence,
)

@pytest.fixture
def simple_system():
    A = np.array(

        [
            [-5,  1,  1,  1],
            [ 1, -5,  1,  1],
            [ 1,  1, -5,  1],
            [ 1,  1,  1, -5]
        ]
    )
    r = np.array([4, 3, 2, 1])
    return A, r


def test_b_is_feasible_Q(simple_system):
    AA, r = simple_system
    assert b_is_feasible_Q(AA, r) == True


def test_simulate_one_step(simple_system):
    nn = 4
    AA = torch.diag(-5*torch.ones(nn))
    rr = 10*torch.ones(nn)
    X0 = torch.arange(nn) + 1
    dt = 1
    X1 = X0 + dt*torch.tensor([5, 0, -15, -40])
    X1 = torch.clamp(X1, min=0) # No negative values allowed
    n_trials = 1
    steps = 1
    sigma = None
    traj = simulate(AA, rr, X0, dt, steps, sigma, n_trials)
    np_traj = traj.cpu().numpy()
    assert traj.shape == (n_trials, steps+1, nn)
    np.testing.assert_allclose(np_traj[0, 0], X0.cpu().numpy(), rtol=1e-5)
    np.testing.assert_allclose(np_traj[0, 1], X1.cpu().numpy(), rtol=1e-5)
    # assert isinstance(traj, torch.Tensor)

def test_random_feasiable_system_x_star():
    nn = 4
    X0 = torch.arange(nn) + 1.0
    # X0 = torch.tensor([0.5, 1.0, 1.5, 2.0])
    device = "cpu"
    AA, rr = random_feasible_system_x_star(X0, device)
    assert rr.shape == X0.shape
    assert AA.shape == (X0.shape[0], X0.shape[0])
    dt = 1
    steps = 10
    sigma = None
    n_trials = 1
    traj = simulate(AA, rr, X0, dt, steps, sigma, n_trials)
    np_traj = traj.cpu().numpy()
    X1 = np_traj[0, 1]
    X10 = np_traj[0, 10]
    np.testing.assert_allclose(np_traj[0, 0], X0.cpu().numpy(), rtol=1e-5)
    np.testing.assert_allclose(X1, X0.cpu().numpy(), rtol=1e-5)
    np.testing.assert_allclose(X10, X0.cpu().numpy(), rtol=1e-5)

@pytest.fixture
def n_species_2_sytem():
    n_species = 2
    X0 = torch.arange(n_species) + 1.0
    # X0 = torch.tensor([0.5, 1.0, 1.5, 2.0])
    device = "cpu"
    A, r = random_feasible_system_x_star(X0, device)
    return n_species, A, r, X0

@pytest.fixture
def n_species_no_GLV():
    n_species = 2
    device = "cpu"
    X0 = torch.arange(n_species, device=device) + 1.0
    A = torch.zeros(n_species, n_species, device=device)
    r = torch.zeros(n_species, device=device)
    return n_species, A, r, X0


@pytest.fixture
def integration_vars():
    dt = 0.1
    steps = 10
    n_trials = 1
    return dt, steps, n_trials

@pytest.fixture
def f_t_n_2_cond_01():
    device = "cpu"
    n_species = 2
    steps = 10
    dt = 0.1
    times = dt*torch.arange(steps, device=device)
    F = torch.zeros(steps, n_species)
    F[:, 0] = 0.1 * torch.sin(times)
    # F[:, 1] = 0.05 * torch.cos(2*times)
    return F, times

@pytest.fixture
def f_t_n_2_cond_02():
    # species 0 with constant perturbation for psteps
    device = "cpu"
    n_species = 2
    steps = 10
    psteps = 5
    dt = 0.1
    times = dt*torch.arange(steps, device=device)
    i_species = 0
    F = torch.zeros((steps, n_species), device=device)
    # f0 = torch.ones(psteps, device=device)
    f0 = 2.0
    F[0:psteps, i_species] = f0
    # F[:, 1] = 0.05 * torch.cos(2*times)
    return F, times


class TestExternalForce:

    def test_simulate_external_force_f_t(self):
        nn = 4
        X0 = torch.arange(nn) + 1.0
        # X0 = torch.tensor([0.5, 1.0, 1.5, 2.0])
        device = "cpu"
        AA, rr = random_feasible_system_x_star(X0, device)
        assert rr.shape == X0.shape
        assert AA.shape == (X0.shape[0], X0.shape[0])
        dt = 1
        steps = 10
        sigma = None
        n_trials = 1
        times = dt*torch.arange(steps, device=device)
        F = torch.zeros(steps, nn)
        F[:, 0] = 0.1 * torch.sin(times)
        F[:, 1] = 0.05 * torch.cos(2 * times)
        # f_t = torch.sin(0.1*time)
        assert F.shape[0] == steps
        traj = simulate(AA, rr, X0, dt, steps, sigma, n_trials, F)
        assert traj.shape == (n_trials, steps+1, nn)

    def test_f_t_wrong_shape(self, n_species_2_sytem, integration_vars, f_t_n_2_cond_01):
        n_species, A, r, X0 = n_species_2_sytem
        dt, steps, n_trials = integration_vars
        F, times = f_t_n_2_cond_01
        assert F.shape[0] == steps
        assert F.shape[1] == n_species
    
    def test_just_f_t_no_GLV(self, n_species_no_GLV, integration_vars, f_t_n_2_cond_01):
        n_species, A, r, X0 = n_species_no_GLV
        dt, steps, n_trials = integration_vars
        F, times = f_t_n_2_cond_01
        sigma = None
        traj = simulate(A, r, X0, dt, steps, sigma, n_trials, F)
        np_traj = traj.cpu().numpy()
        X1 = X0 + dt*F[0]
        np_X0 = X0.cpu().numpy()
        np_X1 = X1.cpu().numpy()
        assert traj.shape == (n_trials, steps+1, n_species)
        np.testing.assert_allclose(np_X0, np_traj[0, 0], rtol=1e-5)
        np.testing.assert_allclose(np_X1, np_traj[0, 1], rtol=1e-5)
    
    def test_f_t_raise_ValueError(self, n_species_no_GLV, integration_vars):
        n_species, A, r, X0 = n_species_no_GLV
        dt, steps, n_trials = integration_vars
        F = torch.zeros(7, 5)
        sigma=None
        with pytest.raises(ValueError):
            traj = simulate(A, r, X0, dt, steps, sigma, n_trials, F)
    
    def test_f_t_simulate(self, n_species_no_GLV, integration_vars, f_t_n_2_cond_02):
        n_species, A, r, X0 = n_species_no_GLV
        dt, steps, n_trials = integration_vars
        F, times = f_t_n_2_cond_02
        assert F.shape == (steps, n_species)