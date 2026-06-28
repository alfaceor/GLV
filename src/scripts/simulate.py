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
    random_feasible_system_x_star,
    random_system,
    simulate
)

from theomodels.plotting import plot_trajectories, plot_matrix_A
from theomodels.io import *
from theomodels.config import Config, register_configs


from theomodels.tracking import start_run, log_artifact_path, end_run

register_configs()  # call before @hydra.main

@hydra.main(version_base=None, config_name="base_config")
def main(cfg: Config) -> None:
    run = start_run(cfg) # Begin tracking experiment
    try:
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
        X0 = torch.arange(cfg.n_species, device=device, dtype=torch.float32) + 1.0
        A, r = random_feasible_system_x_star(X0, device=device)
        # A, r, X0 = random_system(
        #     cfg.n_species,
        #     device=device,
        #     mode=cfg.system_mode,
        # )
        if cfg.save_results:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")


        metadata = {
            "device": str(device),
            "timestamp": timestamp,
            "system_mode": cfg.system_mode,
        }

        output_dir = Path(cfg.output_dir)
        # fln_system = output_dir / f"{cfg.run_name}_n_{cfg.n_species}_system.h5"
        # save_system_to_hdf5(
        #     fln_system,
        #     A.cpu().numpy(),
        #     r.cpu().numpy(),
        #     X0.cpu().numpy(),
        #     )
        
        # 2. External force 
        if cfg.extft.type == "EFNone":
            f_t = None
        else:
            f_t = build_extft(cfg, device=device)
            print(f_t)

        # 3. simulate data
        if cfg.noise.type == "NoneNoise":
            sigma = None
            # traj = simulate(A, r, X0, cfg.dt, cfg.steps)
        else:
            sigma = build_sigma_noise(cfg, device=device)
        print("Simulation started:", datetime.now().strftime("%Y%m%d_%H%M%S"))
        traj = simulate(
            A,
            r,
            X0,
            cfg.dt,
            cfg.steps,
            sigma=sigma,
            n_trials=cfg.n_trials,
            f_t=f_t
        )
        tt = torch.arange(cfg.steps) * cfg.dt
        print("Simulation ended:", datetime.now().strftime("%Y%m%d_%H%M%S"))

        if cfg.save_results:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_dir = Path(cfg.output_dir)
            if cfg.noise.type == "NoneNoise":
                filename = output_dir / f"{cfg.run_name}_n_{cfg.n_species}_traj.h5"
            else:
                filename = output_dir / f"{cfg.run_name}_n_{cfg.n_species}_traj_{cfg.noise.type}.h5"
            metadata = {
                "device": str(device),
                "timestamp": timestamp,
                "system_mode": str(cfg.system_mode),
                "noise": str(cfg.noise)
            }
            sigma_val = sigma.cpu().numpy() if sigma is not None else None
            print(cfg)
            print(type(cfg))

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
    finally:
        end_run() # End tracking experiment



if __name__ == "__main__":
    main()
