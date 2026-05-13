import pytest
import numpy as np
import torch
from theomodels.core import (
    random_feasible_system_x_star,
    b_is_feasible_Q,
    simulate
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

