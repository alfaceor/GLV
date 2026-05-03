"""Hydra entry point that runs the GLV helpers defined in :mod:`theomodels.core`."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import hydra
# from hydra.core.config_store import ConfigStore
from omegaconf import OmegaConf

import torch
import matplotlib.pyplot as plt

from theomodels.core import (
    random_system,
    simulate
)

from theomodels.plotting import plot_trajectories, plot_matrix_A
from theomodels.io import *
from theomodels.config import Config, register_configs


register_configs()  # call before @hydra.main

@hydra.main(version_base=None, config_name="base_config")
def main(cfg: Config) -> None:
    print(OmegaConf.to_yaml(cfg))
    if cfg.seed is not None:
        torch.manual_seed(cfg.seed)

    if cfg.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(cfg.device)

    print("Configuration:", cfg)
    print(f"Using device: {device}")

    # 1. Build system
    A, r, X0 = random_system(
        cfg.n_species,
        device=device,
        mode=cfg.system_mode,
    )
    if cfg.save_results:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")


    metadata = {
        "device": str(device),
        "timestamp": timestamp,
        "system_mode": cfg.system_mode,
    }

    output_dir = Path(cfg.output_dir)
    fln_system = output_dir / f"{cfg.run_name}_n_{cfg.n_species}_system.h5"

    # TODO: Resolve the accordance with the refactored function
    # there should an option to use an already generated file
    # the output should be a file that will be use for a 
    save_system_to_hdf5(
        fln_system,
        A.cpu().numpy(),
        r.cpu().numpy(),
        X0.cpu().numpy(),
        )

    # 2. simulate data
    if cfg.noise._target_ == "NoneNoise":
        sigma = None
        traj = simulate(A, r, X0, cfg.dt, cfg.steps)
    else:
        sigma = build_sigma_noise(cfg, device=device)
        traj = simulate(
            A,
            r,
            X0,
            cfg.dt,
            cfg.steps,
            sigma=sigma,
            n_trials=cfg.n_trials,
        )
    tt = torch.arange(cfg.steps) * cfg.dt

    if cfg.save_results:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path(cfg.output_dir)
        if cfg.noise.type == "NoneNoise":
            filename = output_dir / f"{cfg.run_name}_n_{cfg.n_species}_traj.h5"
        else:
            filename = output_dir / f"{cfg.run_name}_n_{cfg.n_species}_traj_{cfg.noise._target_}.h5"
        metadata = {
            "device": str(device),
            "timestamp": timestamp,
            "system_mode": cfg.system_mode,
            "noise": cfg.noise
        }
        sigma_val = sigma.cpu().numpy() if sigma is not None else None
        
        save_traj_to_hdf5(
            filename,
            traj=traj.cpu().numpy(),
            time=tt.cpu().numpy(),
            A=A.cpu().numpy(),
            r=r.cpu().numpy(),
            X0=X0.cpu().numpy(),
            sigma=sigma_val,
            metadata=metadata,
            cfg=cfg,
        )
        print(f"Saved simulation results to {filename}")



    # fig, ax = plt.subplots(3, 1, figsize=(10, 8))
    # ax_traj, ax_matrix, ax_species_ij = ax

    # for trial in range(cfg.n_trials):
    #     ax_traj.plot(
    #         tt.cpu().numpy(),
    #         stoch_traj[trial, :, 3].cpu().numpy(),
    #         color="gray",
    #         alpha=0.5,
    #         label=f"3, {trial}",
    #     )
    #     ax_species_ij.plot(
    #         stoch_traj[trial, :, 3].cpu().numpy(),
    #         stoch_traj[trial, :, 0].cpu().numpy(),
    #         linestyle="none",
    #         marker="o",
    #     )

    # ax_traj.plot(tt.cpu().numpy(), deter_traj[0, :, 3].cpu().numpy())
    # ax_traj.plot(
    #     deter_traj[0, :, 3].cpu().numpy(),
    #     deter_traj[0, :, 0].cpu().numpy(),
    # )

    # plot_matrix_A(A, ax=ax_matrix)

    # print("Final state (deterministic):", deter_traj[-1].cpu())
    # print("Final state (stochastic):", stoch_traj[-1].cpu())

    # plt.show()


if __name__ == "__main__":
    main()
