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

from theomodels.diversity import compute_diversity, DiversityResults

from theomodels.tracking import start_run, log_artifact_path, end_run

import mlflow

register_configs()  # call before @hydra.main

# ---------------------------------------------------------------------------
# Helper: persist diversity results to the same HDF5 file
# ---------------------------------------------------------------------------

def _save_diversity_to_hdf5(filename: Path, div: DiversityResults) -> None:
    """Append a /diversity group to an existing HDF5 trajectory file.

    Layout
    ------
    /diversity/per_realization/rel_abundance   (n_trials, time, n_species)
    /diversity/per_realization/shannon         (n_trials, time)
    /diversity/per_realization/jsd             (n_trials, time)

    /diversity/ensemble/rel_abundance/{mean,std,q05,q25,median,q75,q95}
    /diversity/ensemble/shannon/{mean,std,q05,q25,median,q75,q95}
    /diversity/ensemble/jsd/{mean,std,q05,q25,median,q75,q95}
    """
    import h5py

    def _write(group, tensor: torch.Tensor) -> None:
        group.create_dataset(
            "data",
            data=tensor.cpu().numpy(),
            compression="gzip",
            compression_opts=4,
        )

    with h5py.File(filename, "a") as f:
        div_grp = f.require_group("diversity")

        # ---- per-realization ----
        pr = div_grp.require_group("per_realization")
        pr.create_dataset("rel_abundance", data=div.rel_abundance.cpu().numpy(),
                          compression="gzip", compression_opts=4)
        pr.create_dataset("shannon", data=div.shannon.cpu().numpy(),
                          compression="gzip", compression_opts=4)
        pr.create_dataset("jsd", data=div.jsd.cpu().numpy(),
                          compression="gzip", compression_opts=4)

        # ---- ensemble statistics ----
        ens = div_grp.require_group("ensemble")
        for metric_name, stats in (
            ("rel_abundance", div.rel_abundance_stats),
            ("shannon",       div.shannon_stats),
            ("jsd",           div.jsd_stats),
        ):
            grp = ens.require_group(metric_name)
            for stat_name, tensor in stats.to_dict().items():
                grp.create_dataset(
                    stat_name,
                    data=tensor.cpu().numpy(),
                    compression="gzip",
                    compression_opts=4,
                )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

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

        # ------------------------------------------------------------------
        # 1. Build system
        # ------------------------------------------------------------------
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
        
        # ------------------------------------------------------------------
        # 2. External force
        # ------------------------------------------------------------------
        if cfg.extft.type == "EFNone":
            f_t = None
        else:
            f_t = build_extft(cfg, device=device)

        # ------------------------------------------------------------------
        # 3. Simulate   traj: (n_trials, steps+1, n_species)
        # ------------------------------------------------------------------
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
        tt = torch.arange(cfg.steps+1) * cfg.dt
        print("Simulation ended:", datetime.now().strftime("%Y%m%d_%H%M%S"))

        # ------------------------------------------------------------------
        # 4. Diversity metrics
        #    - relative abundance : (n_trials, steps+1, n_species)
        #    - Shannon index      : (n_trials, steps+1)
        #    - JSD vs ensemble    : (n_trials, steps+1)
        #    Each carries ensemble mean / std / quantiles over trials.
        # ------------------------------------------------------------------
        print("Computing diversity metrics...")
        div = compute_diversity(traj, base=None)   # nats; use base=2 for bits

        print(
            f"  Shannon — mean at t=0: {div.shannon_stats.mean[0]:.4f} nats, "
            f"mean at t=final: {div.shannon_stats.mean[-1]:.4f} nats"
        )
        print(
            f"  JSD vs ensemble — mean at t=0: {div.jsd_stats.mean[0]:.4f}, "
            f"mean at t=final: {div.jsd_stats.mean[-1]:.4f}"
        )

                # ── log anything you like, directly, no boilerplate needed ──
        mlflow.log_metric("diversity.shannon.final_mean",   div.shannon_stats.mean[-1].item())
        mlflow.log_metric("diversity.shannon.final_std",    div.shannon_stats.std[-1].item())
        mlflow.log_metric("diversity.shannon.final_q25",    div.shannon_stats.q25[-1].item())
        mlflow.log_metric("diversity.shannon.final_median", div.shannon_stats.median[-1].item())
        mlflow.log_metric("diversity.shannon.final_q75",    div.shannon_stats.q75[-1].item())

        mlflow.log_metric("diversity.jsd.final_mean",       div.jsd_stats.mean[-1].item())
        mlflow.log_metric("diversity.jsd.final_std",        div.jsd_stats.std[-1].item())


        # ------------------------------------------------------------------
        # 5. Persist
        # ------------------------------------------------------------------
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
