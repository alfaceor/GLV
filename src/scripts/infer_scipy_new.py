# %%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from pathlib import Path

from theomodels.io import load_extft_from_hdf5, load_traj_from_hdf5
from theomodels.paths import DATA
from theomodels.infer_glv import GLVInferenceEngine
# Ensure plots are clean and readable
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')


# Execution Block
if __name__ == "__main__":

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


    
    # # Mocking standard inputs to replace your file loader logic for isolated validation testing
    # time_mock = np.linspace(0, 10, 50)
    # # Generate clean toy trajectory: N=2 species
    # true_r = np.array([5.0, 3.0])
    # true_A = np.array([[-4.0, -1.5], [-1.0, -3.0]])
    
    # def toy_system(t, X):
    #     return X * (true_r + true_A @ X)
        
    # sol_toy = solve_ivp(toy_system, (0, 10), y0=[1.0, 0.8], t_eval=time_mock)
    # observed_mock = sol_toy.y.T + np.random.normal(0, 0.05, sol_toy.y.T.shape)

    # Instantiate Inference Engine
    engine = GLVInferenceEngine(time_data=time2infer, trajectory_data=data2infer)

    ini_r = np.array([8, 3.5])
    # ini_A = -1.0*np.eye(nn).flatten()
    # ini_A = np.array([-4., -2., -2., -1.])
    ini_A = np.array([-4.8, -1.8, -1.8, -0.9])
    # ini_A = np.array([-4., -2., -1, -1.])
    # ini_r = data_r
    # ini_A = data_A.flatten()

    ini_sigma = np.array([1.0])
    param_initial_guess = np.concatenate([ini_r, ini_A])
    param_initial_guess = np.concatenate([param_initial_guess, ini_sigma])

    
    # # Define perturbation parameters manually if desired
    # manual_initial_guess = np.concatenate([
    #     [4.0, 2.5],             # r guess
    #     [-3.0, -1.0, -1.0, -2.0], # A matrix guess elements flat
    #     [1.0]                   # Sigma noise guess
    # ])

    # Run fitting
    fit_res = engine.fit(custom_initial_guesses=param_initial_guess)
    print("\n--- INFERENCE RESULTS ---")
    print(f"Convergence State: {fit_res.success}")
    print(f"Estimated parameters:\n {fit_res.x}")

    r_inferred = fit_res.x[:engine.num_species]
    A_inferred = fit_res.x[engine.num_species:engine.num_species*(engine.num_species+1)].reshape(engine.num_species, engine.num_species)
    
    inferred_traj = engine.integrate(time2infer, r_inferred, A_inferred, y0=None)

    np_engine_nll_history = np.array(engine.nll_history)
    np_engine_nll_history[np_engine_nll_history > 1e4] = np.nan
    df_history = pd.DataFrame(engine.param_history, columns=['r1', 'r2', 'a11', 'a12', 'a21', 'a22', 'sigma'])
    df_history["nll"] = engine.nll_history
    df_history

    # Plot optimization diagnostics
    fig, axes = plt.subplots(2, 2, figsize=(12, 5))
    axnll = axes[0][0]
    axdata = axes[1][0]
    axparall = axes[0][1]

    axnll.plot(np_engine_nll_history, 'o', color='crimson', lw=2)
    axnll.set_title("Negative Log-Likelihood Improvement")
    axnll.set_xlabel("Function Evaluation Steps")
    axnll.set_ylabel("NLL")
    # axnll.set_yscale("log")

    axdata.plot(inferred_traj)
    axdata.plot(data2infer, 'o')

    pd.plotting.parallel_coordinates(df_history, ax=axparall, class_column='nll', colormap='Reds', linewidth=2, legend=False)

    # Generate Pairplot of Parameter Evolution trajectories
    
    # Take every 5th iteration step to prevent rendering lag
    # sns.pairplot(df_history.iloc[::5, :], diag_kind='kde', corner=True)
    # plt.show()
    import seaborn as sns
    
    # import plotly.express as px
    # px.parallel_coordinates(df_history, color='nll')
    plt.show()

# %%
