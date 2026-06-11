# %%
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
font = {'family': 'normal',
        'weight': 'bold',
        'size': 14}
plt.rc('font', **font)
plt.rc('text', usetex=True)

from theomodels.io import load_extft_from_hdf5, load_traj_from_hdf5
from theomodels.paths import DATA

from scipy.integrate import solve_ivp
from scipy.optimize import minimize

from sklearn.preprocessing import StandardScaler



# 1. Load data
PATH_DATA_SIMUL = DATA / "simul"

PATH_ALL_EFS = PATH_DATA_SIMUL / "glv_all_efs_n_2_traj.h5"
PATH_0_EFS_NA = PATH_DATA_SIMUL / "glv_0_efs_na_n_2_traj_AllNoise.h5"
PATH_1_EFS_NA = PATH_DATA_SIMUL / "glv_1_efs_na_n_2_traj_AllNoise.h5"
PATH_ALL_EFS = PATH_DATA_SIMUL / "glv_all_efs_n_2_traj.h5"


data_0_efs_na = load_traj_from_hdf5(PATH_0_EFS_NA)
data_1_efs_na = load_traj_from_hdf5(PATH_1_EFS_NA)


# selected data
data = data_1_efs_na
extft_index = 1
i_trial = 10

n_trials, n_time, n_species = data["traj"].shape

data_A = data['system']['A']
data_r = data['system']['r']
data_X0 = data['system']['X0']

data_traj = data["traj"]
data_sys = data["traj"]
time = data["time"]


# 2. Select samples
aaa = 100
data2infer = data_traj[i_trial, :-1:aaa, :]
time2infer = time[::aaa]
data2infer.shape, time2infer.shape


# 3. Define ode system
def ode_system(t, X, r, A, tpert=0.2, f_0=-4.0):
    # f_0 = -4.0
    f_t = np.zeros_like(r)
    if t < tpert:
        f_t[extft_index] = f_0
    dydt = X * (r + A @ X) + f_t
    return dydt


# 4. Define negative log likelihood
def negative_log_likelihood(params_to_estimate, t_eval, observed_data, y0):
    # Unpack parameters we are trying to fit
    da_r = params_to_estimate[0:2]
    da_A = params_to_estimate[2:6].reshape(2, 2)
    sigma = params_to_estimate[6]

    # Solve the ODE with the current parameter guesses
    sol = solve_ivp(
        ode_system,
        t_span=(t_eval[0], t_eval[-1]),
        y0=y0,
        t_eval=t_eval,
        args=(da_r, da_A)
    )
    # If the ODE solver failed to converge, return a high penalty
    if not sol.success:
        return np.inf

    # Calculate the Negative Log-Likelihood for Gaussian noise
    # Formula: N * log(sigma) + sum((obs - pred)^2) / (2 * sigma^2)
    n = observed_data.size
    residuals = np.log(observed_data) - np.log(sol.y.T)

    # FIXME: Change likelihood
    nll = (n / 2) * np.log(2 * np.pi) + n * np.log(sigma) + np.sum(residuals**2) / (2 * sigma**2)
    return nll



# 5. Create an empty list to store the history
nll_history = []
xk_history = []

# 6. Define the callback function
# SciPy automatically passes the current parameter array (xk) to this function
def track_evolution(xk):
    # Calculate the NLL for the current parameters and store it
    # We reuse your existing negative_log_likelihood function
    current_nll = negative_log_likelihood(xk, time2infer, data2infer, data_X0)
    nll_history.append(current_nll)
    xk_history.append(xk)


# Define your custom limits
custom_options = {
    'maxiter': 2000,  # Max allowed iterations (Default is usually N * 200)
    'maxfev': 5000,   # Max allowed function evaluations
    'disp': True      # Set to True to print convergence messages right in the console
}

# 7. Minimize log-likelihood
nn = 2
# FIXME: The initial conditions of r play a fundamental role to obtain a reasonable fit of the data
# TODO: I would like to make a perturbation of the parameters and explore the influence and the uncertainty of the fit. Then how can I get a better guess of the initial conditions what should be the tip & tricks for that.
ini_r = np.array([8, 3.5])
# ini_A = -1.0*np.eye(nn).flatten()
# ini_A = np.array([-4., -2., -2., -1.])
# ini_A = np.array([-4.8, -1.8, -1.8, -0.9])
ini_A = np.array([-4., -2., -1, -1.])
# ini_r = data_r
# ini_A = data_A.flatten()

ini_sigma = np.array([1.0])
initial_guess = np.concatenate([ini_r, ini_A])
initial_guess = np.concatenate([initial_guess, ini_sigma])
# initial_guess = [ini_r, ini_A]

# Run the optimization
result = minimize(
    negative_log_likelihood, 
    initial_guess, 
    args=(time2infer, data2infer, data_X0),
    method='Nelder-Mead', # 'Nelder-Mead' or 'L-BFGS-B' work well for ODEs
    callback=track_evolution,
    # options=custom_options
)

# Extract fitted parameters
fitted_params = result.x

rr = fitted_params[0:2]
AA = fitted_params[2:6].reshape(2, 2)
ss = fitted_params[6]

sol_infer = solve_ivp(
    ode_system,
    t_span=(time[0], time[-1]),
    y0=data_X0,
    t_eval=time,
    args=(rr, AA)
)

sol = solve_ivp(
        ode_system, 
        t_span=(time[0], time[-1]), 
        y0=data_X0, 
        t_eval=time, 
        args=(data_r, data_A)
    )

fig, ax = plt.subplots(2, 2)

ax_fits = ax[0][0]
ax_nll = ax[0][1]
ax_xk = ax[1][1]
ax_sca = ax[1][0]

ax_fits.plot(sol.t, sol.y.T, '--')
ax_fits.plot(time2infer, data2infer, 'o')
ax_fits.plot(sol_infer.t, sol_infer.y.T)


ax_nll.plot(nll_history)
np_xk_history = np.array(xk_history)
print(np_xk_history.shape)

import seaborn as sns
import pandas as pd
df = pd.DataFrame(
    np_xk_history,
    columns=['r_1', 'r_2', 'a00', 'a01', 'a10', 'a11', 'sigma']
)
sns.pairplot(df)


theo_params = np.concatenate(
    [
        data_r,
        data_A.flatten()
    ]
)

theo_params = np.concatenate(
    [
        theo_params,
        np.array([ data["config"]["noise"]["std"] ])
    ]
)

ax_xk.plot(
    xk_history, 'o'
)

print('-'*50)
print(theo_params)
print(fitted_params)
print('-'*50)

data["config"]


from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
np_xk_history_scaler = scaler.fit_transform(np_xk_history)

ax_sca.plot(np_xk_history_scaler)

import plotly.express as px


# px.parallel_coordinates()

plt.show()

# print("Optimization Successful:", result.success)
# print(f"Fitted alpha, beta, delta, gamma: {fitted_params[:4]}")
# print(f"Estimated noise (sigma): {fitted_params[4]:.4f}")




