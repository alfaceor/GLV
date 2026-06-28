"""Query and load simulation runs from MLflow tracking."""
from __future__ import annotations
import mlflow
import h5py
import numpy as np
import pandas as pd
from theomodels.tracking import TRACKING_URI, EXPERIMENT


def get_client() -> mlflow.MlflowClient:
    mlflow.set_tracking_uri(TRACKING_URI)
    return mlflow.MlflowClient()


def query_runs(filters: dict[str, str | float | list]) -> pd.DataFrame:
    """
    Query runs by parameter values.

    filters keys match flattened param names logged in tracking.py.
    Pass a list for OR matching (e.g. std in [0.1, 0.5, 1.0]).

    Example
    -------
    query_runs({
        "noise.type":  "SingleNoise",
        "extft.type":  "EFSingle",
        "extft.f0":    -2.0,
        "noise.index": 1,
        "noise.std":   [0.1, 0.5, 1.0],
    })
    """
    client = get_client()
    exp = client.get_experiment_by_name(EXPERIMENT)
    if exp is None:
        raise ValueError(f"Experiment '{EXPERIMENT}' not found.")

    # Build MLflow filter string for scalar params
    clauses = []
    list_filters = {}
    for k, v in filters.items():
        if isinstance(v, list):
            list_filters[k] = [str(x) for x in v]
        else:
            clauses.append(f'params.`{k}` = "{v}"')

    filter_str = " and ".join(clauses)
    runs = client.search_runs(
        experiment_ids=[exp.experiment_id],
        filter_string=filter_str,
        max_results=500,
    )

    # Post-filter list params (MLflow doesn't support IN natively)
    if list_filters:
        runs = [
            r for r in runs
            if all(
                r.data.params.get(k) in vals
                for k, vals in list_filters.items()
            )
        ]

    rows = []
    for r in runs:
        row = {**r.data.params, **r.data.tags}
        row["run_id"]    = r.info.run_id
        row["status"]    = r.info.status
        row["start_time"] = pd.Timestamp(r.info.start_time, unit="ms")
        rows.append(row)

    return pd.DataFrame(rows)


def load_traj(run_row: pd.Series) -> dict:
    """Load HDF5 arrays for a single run from a query_runs() result row."""
    path = run_row["h5_path"]
    with h5py.File(path, "r") as f:
        return {
            "traj":  f["traj"][:],    # (n_trials, steps, n_species)
            "time":  f["time"][:],
            "A":     f["A"][:],
            "r":     f["r"][:],
            "X0":    f["X0"][:],
            "sigma": f["sigma"][:] if "sigma" in f else None,
        }


def load_runs(filters: dict) -> tuple[pd.DataFrame, list[dict]]:
    """Convenience: query + load arrays in one call.
    
    Returns (metadata_df, list_of_data_dicts) in the same row order.
    """
    df = query_runs(filters)
    data = [load_traj(row) for _, row in df.iterrows()]
    return df, data