"""Thin MLflow wrapper for GLV simulation tracking."""
from __future__ import annotations
import mlflow
from omegaconf import OmegaConf
from pathlib import Path

TRACKING_URI = "sqlite:///mlflow.db"   # rsync this file from cluster
# TRACKING_URI = "sqlite:///data/mlruns.db"   # rsync this file from cluster
EXPERIMENT   = "glv_sims"  # Generalized lotka-volterra simulations


def start_run(cfg) -> mlflow.ActiveRun:
    mlflow.set_tracking_uri(TRACKING_URI)
    mlflow.set_experiment(EXPERIMENT)
    params = OmegaConf.to_container(cfg, resolve=True, throw_on_missing=False)
    # Flatten nested dicts (noise.std, extft.f0, …)
    flat = _flatten(params)
    run = mlflow.start_run(run_name=cfg.run_name)
    mlflow.log_params(flat)
    return run


def log_artifact_path(h5_path: str | Path) -> None:
    """Store the absolute path as a tag — avoids copying large files."""
    mlflow.set_tag("h5_path", str(Path(h5_path).resolve()))


def end_run() -> None:
    mlflow.end_run()


def _flatten(d: dict, parent: str = "") -> dict:
    out = {}
    for k, v in d.items():
        key = f"{parent}.{k}" if parent else k
        if isinstance(v, dict):
            out.update(_flatten(v, key))
        else:
            out[key] = v
    return out