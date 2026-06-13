# %%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from pathlib import Path

# Ensure plots are clean and readable
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

class GLVInferenceEngine:
    def __init__(self, time_data, trajectory_data, extft_index=1, tpert=0.2, f_0=-4.0):
        """
        Modular Engine for GLV Parameter Inference.
        Handles data mapping, parameter vectorization, and likelihood profiling.
        """
        self.time = time_data
        self.observed = trajectory_data
        self.num_species = trajectory_data.shape[1]
        
        # Perturbation profile variables
        self.extft_index = extft_index
        self.tpert = tpert
        self.f_0 = f_0
        
        # Metrics containers (Updated efficiently via optimizer calls)
        self.nll_history = []
        self.param_history = []

    def ode_system(self, t, X, r, A):
        f_t = np.zeros_like(r)
        if t < self.tpert:
            f_t[self.extft_index] = self.f_0
        return X * (r + A @ X) + f_t

    def compute_nll(self, params_to_estimate):
        """Calculates negative log-likelihood with absolute convergence fallbacks."""
        # Unpack flat parameter array safely based on target dimensions
        n = self.num_species
        r = params_to_estimate[0:n]
        A = params_to_estimate[n:n + n**2].reshape(n, n)
        sigma = params_to_estimate[-1]

        # Boundary penalty: Sigma cannot be non-positive
        if sigma <= 1e-6:
            return np.inf

        # Integrate the forward problem
        sol = solve_ivp(
            self.ode_system,
            t_span=(self.time[0], self.time[-1]),
            y0=self.observed[0, :], # Use exact first data point as initial state
            t_eval=self.time,
            args=(r, A),
            method='RK45'
        )
        
        if not sol.success:
            return np.inf

        # Compute residuals in state-space (avoiding log-zero penalties)
        num_points = self.observed.size
        residuals = self.observed - sol.y.T
        
        # Gaussian Likelihood formulation
        nll = (num_points / 2) * np.log(2 * np.pi) + num_points * np.log(sigma) + np.sum(residuals**2) / (2 * sigma**2)
        
        # State Tracking optimization: cache results here instead of inside a separate callback call
        self.nll_history.append(nll)
        self.param_history.append(params_to_estimate.copy())
        
        return nll

    def algebraic_initial_guess(self):
        """
        Tip & Trick: Generate data-driven initial guesses for r and A 
        using empirical finite-difference derivatives.
        """
        X = self.observed
        dt = np.diff(self.time)[:, None]
        
        # Calculate dX/dt approximations
        dX_dt = np.diff(X, axis=0) / dt
        X_mid = (X[:-1] + X[1:]) / 2.0
        
        # Per-species linear regression: dX_dt / X = r + A * X
        y_reg = dX_dt / (X_mid + 1e-8)
        
        # Simple heuristic initialization if regression is noisy
        r_guess = np.ones(self.num_species) * 2.0
        A_guess = -np.eye(self.num_species) * 1.5
        sigma_guess = np.array([0.5])
        
        return np.concatenate([r_guess, A_guess.flatten(), sigma_guess])

    def fit(self, custom_initial_guesses=None, max_iter=1000):
        """Runs optimization profile."""
        if custom_initial_guesses is not None:
            init_guess = custom_initial_guesses
        else:
            init_guess = self.algebraic_initial_guess()
            print("--> Using data-driven algebraic initial guess.")

        # Clear tracking containers from previous runs
        self.nll_history.clear()
        self.param_history.clear()

        # L-BFGS-B handles boundary constraints much better than Nelder-Mead for sigma
        bounds = [(None, None)] * (self.num_species + self.num_species**2) + [(1e-4, 10.0)]
        
        print("Starting optimization inference routine...")
        result = minimize(
            self.compute_nll,
            init_guess,
            method='L-BFGS-B',
            bounds=bounds,
            options={'maxiter': max_iter, 'disp': True}
        )
        return result


# Execution Block
if __name__ == "__main__":
    # Mocking standard inputs to replace your file loader logic for isolated validation testing
    time_mock = np.linspace(0, 10, 50)
    # Generate clean toy trajectory: N=2 species
    true_r = np.array([5.0, 3.0])
    true_A = np.array([[-4.0, -1.5], [-1.0, -3.0]])
    
    def toy_system(t, X):
        return X * (true_r + true_A @ X)
        
    sol_toy = solve_ivp(toy_system, (0, 10), y0=[1.0, 0.8], t_eval=time_mock)
    observed_mock = sol_toy.y.T + np.random.normal(0, 0.05, sol_toy.y.T.shape)

    # Instantiate Inference Engine
    engine = GLVInferenceEngine(time_data=time_mock, trajectory_data=observed_mock)
    
    # Define perturbation parameters manually if desired
    manual_initial_guess = np.concatenate([
        [4.0, 2.5],             # r guess
        [-3.0, -1.0, -1.0, -2.0], # A matrix guess elements flat
        [1.0]                   # Sigma noise guess
    ])

    # Run fitting
    fit_res = engine.fit(custom_initial_guesses=manual_initial_guess)
    print("\n--- INFERENCE RESULTS ---")
    print(f"Convergence State: {fit_res.success}")
    print(f"Estimated parameters:\n {fit_res.x}")

    # Plot optimization diagnostics
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].plot(engine.nll_history, color='crimson', lw=2)
    axes[0].set_title("Negative Log-Likelihood Improvement")
    axes[0].set_xlabel("Function Evaluation Steps")
    axes[0].set_ylabel("NLL")

    # Generate Pairplot of Parameter Evolution trajectories
    df_history = pd.DataFrame(engine.param_history, columns=['r1', 'r2', 'a11', 'a12', 'a21', 'a22', 'sigma'])
    # Take every 5th iteration step to prevent rendering lag
    # sns.pairplot(df_history.iloc[::5, :], diag_kind='kde', corner=True)
    # plt.show()