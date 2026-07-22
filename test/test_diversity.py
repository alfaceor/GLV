import torch
import pytest
from theomodels.diversity import (
    relative_abundance,
    shannon_diversity,
    kl_divergence,
    js_divergence,
)


# ---------------------------------------------------------------------------
# relative_abundance
# ---------------------------------------------------------------------------

class TestRelativeAbundance:
    def test_sums_to_one(self):
        X = torch.tensor([1.0, 2.0, 3.0, 4.0])
        p = relative_abundance(X)
        assert torch.allclose(p.sum(), torch.tensor(1.0), atol=1e-6)

    def test_uniform_input(self):
        X = torch.tensor([5.0, 5.0, 5.0, 5.0])
        p = relative_abundance(X)
        expected = torch.full((4,), 0.25)
        assert torch.allclose(p, expected, atol=1e-6)

    def test_known_values(self):
        X = torch.tensor([1.0, 1.0, 2.0])
        p = relative_abundance(X)
        expected = torch.tensor([0.25, 0.25, 0.5])
        assert torch.allclose(p, expected, atol=1e-6)

    def test_batched_input(self):
        X = torch.tensor([[1.0, 1.0, 2.0], [2.0, 2.0, 4.0]])
        p = relative_abundance(X)
        assert p.shape == X.shape
        sums = p.sum(dim=-1)
        assert torch.allclose(sums, torch.ones(2), atol=1e-6)

    def test_negative_values_raise(self):
        X = torch.tensor([1.0, -1.0, 2.0])
        with pytest.raises(ValueError):
            relative_abundance(X)

    def test_zero_vector_does_not_nan(self):
        X = torch.zeros(3)
        p = relative_abundance(X)
        assert torch.all(torch.isfinite(p))

    def test_single_species(self):
        X = torch.tensor([7.0])
        p = relative_abundance(X)
        assert torch.allclose(p, torch.tensor([1.0]), atol=1e-6)


# ---------------------------------------------------------------------------
# shannon_diversity
# ---------------------------------------------------------------------------

class TestShannonDiversity:
    def test_uniform_distribution_is_maximal(self):
        n = 4
        X = torch.ones(n)
        H = shannon_diversity(X)
        assert torch.allclose(H, torch.log(torch.tensor(n)), atol=1e-4)

    def test_single_species_is_zero(self):
        X = torch.tensor([1.0, 0.0, 0.0])
        H = shannon_diversity(X, eps=1e-12)
        assert H.item() == pytest.approx(0.0, abs=1e-6)

    def test_non_negative(self):
        X = torch.tensor([3.0, 1.0, 0.5, 2.5])
        H = shannon_diversity(X)
        assert H.item() >= 0.0

    def test_base_2_vs_natural_log(self):
        X = torch.tensor([1.0, 1.0, 1.0, 1.0])
        H_nat = shannon_diversity(X)
        H_bits = shannon_diversity(X, base=2)
        assert torch.allclose(H_nat / torch.log(torch.tensor(2)), H_bits, atol=1e-4)

    def test_two_equal_species_natural_log(self):
        X = torch.tensor([1.0, 1.0])
        H = shannon_diversity(X)
        assert H.item() == pytest.approx(torch.log(torch.tensor(2)), abs=1e-5)

    def test_batched_trajectory(self):
        X = torch.rand(5, 10, 6) + 0.1  # (n_trials, steps, n_species)
        H = shannon_diversity(X)
        assert H.shape == (5, 10)
        assert torch.all(H >= 0)


# ---------------------------------------------------------------------------
# kl_divergence
# ---------------------------------------------------------------------------

class TestKLDivergence:
    def test_self_divergence_is_zero(self):
        X = torch.tensor([1.0, 2.0, 3.0, 4.0])
        D = kl_divergence(X, X)
        assert D.item() == pytest.approx(0.0, abs=1e-5)

    def test_non_negative(self):
        P = torch.tensor([1.0, 2.0, 3.0])
        Q = torch.tensor([3.0, 1.0, 1.0])
        D = kl_divergence(P, Q)
        assert D.item() >= -1e-6

    def test_asymmetry(self):
        P = torch.tensor([0.6, 0.3, 0.1])
        Q = torch.tensor([0.2, 0.2, 0.6])
        D_pq = kl_divergence(P, Q, from_abundances=False)
        D_qp = kl_divergence(Q, P, from_abundances=False)
        assert not torch.allclose(D_pq, D_qp, atol=1e-3)

    def test_known_value_two_point_distributions(self):
        # D_KL(P||Q) for P=[0.5,0.5], Q=[0.9,0.1] computed by hand (nats)
        P = torch.tensor([0.5, 0.5])
        Q = torch.tensor([0.9, 0.1])
        expected = 0.5 * torch.log(torch.tensor(0.5 / 0.9)) + 0.5 * torch.log(torch.tensor(0.5 / 0.1))
        D = kl_divergence(P, Q, from_abundances=False, eps=0.0)
        assert D.item() == pytest.approx(expected, abs=1e-4)

    def test_from_abundances_matches_normalized(self):
        P_raw = torch.tensor([2.0, 6.0])
        Q_raw = torch.tensor([5.0, 5.0])
        D_raw = kl_divergence(P_raw, Q_raw, from_abundances=True)

        P_norm = relative_abundance(P_raw)
        Q_norm = relative_abundance(Q_raw)
        D_norm = kl_divergence(P_norm, Q_norm, from_abundances=False)

        assert torch.allclose(D_raw, D_norm, atol=1e-5)

    def test_base_2_scaling(self):
        P = torch.tensor([0.7, 0.3])
        Q = torch.tensor([0.4, 0.6])
        D_nat = kl_divergence(P, Q, from_abundances=False)
        D_bits = kl_divergence(P, Q, from_abundances=False, base=2)
        assert torch.allclose(D_nat / torch.log(torch.tensor(2)), D_bits, atol=1e-4)

    def test_batched_shapes(self):
        P = torch.rand(4, 5) + 0.1
        Q = torch.rand(4, 5) + 0.1
        D = kl_divergence(P, Q)
        assert D.shape == (4,)


# ---------------------------------------------------------------------------
# js_divergence
# ---------------------------------------------------------------------------

class TestJSDivergence:
    def test_self_divergence_is_zero(self):
        X = torch.tensor([1.0, 2.0, 3.0, 4.0])
        D = js_divergence(X, X)
        assert D.item() == pytest.approx(0.0, abs=1e-5)

    def test_symmetry(self):
        P = torch.tensor([0.8, 0.1, 0.1])
        Q = torch.tensor([0.1, 0.1, 0.8])
        D_pq = js_divergence(P, Q, from_abundances=False)
        D_qp = js_divergence(Q, P, from_abundances=False)
        assert torch.allclose(D_pq, D_qp, atol=1e-5)

    def test_non_negative(self):
        P = torch.tensor([1.0, 2.0, 3.0])
        Q = torch.tensor([3.0, 1.0, 1.0])
        D = js_divergence(P, Q)
        assert D.item() >= -1e-6

    def test_upper_bound_natural_log(self):
        # JSD is bounded above by log(2) (nats) for any two distributions
        P = torch.tensor([1.0, 0.0])
        Q = torch.tensor([0.0, 1.0])
        D = js_divergence(P, Q, from_abundances=False, eps=1e-12)
        assert D.item() <= torch.log(torch.tensor(2)) + 1e-4

    def test_disjoint_supports_near_log2(self):
        P = torch.tensor([1.0, 1e-12])
        Q = torch.tensor([1e-12, 1.0])
        D = js_divergence(P, Q, from_abundances=False, eps=1e-12)
        assert D.item() == pytest.approx(torch.log(torch.tensor(2)), abs=1e-2)

    def test_base_2_in_unit_interval(self):
        P = torch.tensor([1.0, 0.0])
        Q = torch.tensor([0.0, 1.0])
        D = js_divergence(P, Q, from_abundances=False, base=2)
        assert D.item() <= 1.0 + 1e-4

    def test_from_abundances_matches_normalized(self):
        P_raw = torch.tensor([2.0, 6.0, 2.0])
        Q_raw = torch.tensor([1.0, 1.0, 8.0])
        D_raw = js_divergence(P_raw, Q_raw, from_abundances=True)

        P_norm = relative_abundance(P_raw)
        Q_norm = relative_abundance(Q_raw)
        D_norm = js_divergence(P_norm, Q_norm, from_abundances=False)

        assert torch.allclose(D_raw, D_norm, atol=1e-5)

    def test_batched_shapes(self):
        P = torch.rand(3, 4, 6) + 0.1
        Q = torch.rand(3, 4, 6) + 0.1
        D = js_divergence(P, Q)
        assert D.shape == (3, 4)

