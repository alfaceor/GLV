import torch
# from dataclasses import dataclass, field
# from typing import Any
from theomodels.config import Config
import h5py
import numpy as np
from pathlib import Path


def resolve_device(device_str: str) -> torch.device:
    if device_str == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(device_str)


def parse_noise_map(noise_map_str: str, n_species: int) -> dict[int, float]:
    """Parse '0:0.1;2:0.4' into {0: 0.1, 2: 0.4} — semicolon separated"""
    if not noise_map_str.strip():
        return {}
    
    noise_map = {}
    for pair in noise_map_str.split(";"):        # ← semicolon separator
        parts = pair.strip().split(":")
        if len(parts) != 2:
            raise ValueError(f"Invalid pair '{pair}', expected format 'key:value'")
        key, value = int(parts[0].strip()), float(parts[1].strip())
        if not (0 <= key < n_species):
            raise ValueError(f"Key {key} out of range [0, {n_species})")
        noise_map[key] = value
    
    return noise_map



def build_sigma_all(cfg: Config, device: torch.device) -> torch.Tensor:
    sigma = torch.full((cfg.n_species,), cfg.noise_std, device=device)
    return sigma

def build_sigma_single(cfg: Config, device: torch.device) -> torch.Tensor:
    sigma = torch.zeros(cfg.n_species, device=device)
    sigma[cfg.index] = cfg.std
    return sigma
    # raise NotImplementedError


def build_sigma_selid(cfg: Config, device: torch.device) -> torch.Tensor:
    # generate dictionary

    sigma = torch.full((cfg.n_species,), cfg.default_std, device=device)
    # TODO: Finish implementation:    
    base_std = cfg.noise_std
    noise_map = parse_noise_map(cfg.noise_map_str, cfg.n_species)
    try:
        for key, value in cfg.noise.noise_map.items():
            sigma[key] = value
    except Exception as e:
        print(e)
        raise e
    return sigma


def build_sigma_noise(cfg: Config, device: torch.device) -> torch.Tensor:
    """Return the sigma tensor with the noise intensity for each species

    Args:
        cfg (Config): Hydra configuration instance of inputs. See theomodels.config.py

    Raises:
        NotImplementedError: _description_

    Returns:
        torch.Tensor: Noise intensity tensor
    """
    if cfg.noise._target_ == "AllNoise":
        sigma = build_sigma_all(cfg, device)
    elif cfg.noise._target_ == "single":
        sigma = build_sigma_selid(cfg, device)
    elif cfg.noise._target_ == "selid":
        sigma = build_sigma_selid(cfg, device)
    elif cfg.noise._target_ == "none":
        sigma = None
    else:
        raise ValueError("Unknow noise configuration")
    return sigma




# def build_sigma(cfg: Config, device: torch.device) -> torch.Tensor:
#     noise_target = cfg.noise._target_

#     if noise_target == "AllNoise":
#         sigma = torch.full((cfg.n_species,), cfg.noise.std, device=device)

#     elif noise_target == "SelectedNoise":
#         sigma = torch.full((cfg.n_species,), cfg.noise.default_std, device=device)
#         noise_map = parse_noise_map(cfg.noise.noise_map_str, cfg.n_species)
#         for key, value in noise_map.items():
#             sigma[key] = value

#     else:
#         raise ValueError(f"Unknown noise target: '{noise_target}'")

#     return sigma




def save_system_to_hdf5(
    filepath: Path,
    cfg: Config,
    A: torch.Tensor,
    r: torch.Tensor,
    X0: torch.Tensor,
    metadata: dict[str, str],
) -> None:
    """ Save growth rate, interaction matrices with implicit carrying capacities

    Args:
        filepath (Path): filepath to save data in h5 file format
        cfg (Config): configuration data
        A (torch.Tensor): Interaction matrix between species
        r (torch.Tensor): growth rate
        X0 (torch.Tensor): initial conditions
        metadata (dict[str, str]): system creation metadata
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(filepath, "w") as h5:
        h5.create_dataset("system/A", data=A.cpu().numpy())
        h5.create_dataset("system/r", data=r.cpu().numpy())
        h5.create_dataset("system/X0", data=X0.cpu().numpy())


def explore(name, obj):
    print(name)
    if isinstance(obj, h5py.Dataset):
        print("  shape:", obj.shape)
        print("  dtype:", obj.dtype)
    print("  attrs:", dict(obj.attrs))


def load_simulation_data(filepath):
    with h5py.File(filepath, "r") as h5:
        # 1. Read Datasets (as NumPy arrays)
        # Use [:] to load the data into memory
        data = {
            "det": h5["trajectories/deterministic"][:],
            "stoch": h5["trajectories/stochastic"][:],
            "A": h5["system/A"][:],
            "r": h5["system/r"][:],
            "X0": h5["system/X0"][:],
            "sigma": h5["system/sigma"][:]
        }

        # 2. Read Config Attributes
        # .attrs returns an object you can cast to a dict
        config = dict(h5["config"].attrs)
        # config = ast.literal_eval(h5["config"].attrs["_content"])

        # 3. Read Metadata Attributes
        metadata = dict(h5["metadata"].attrs)

        return data, config, metadata


def load_system_from_hdf5(filepath: Path) -> dict[str, np.ndarray]:
    """Load system data from hdf5

    Args:
        filepath (Path): filepath of hdf5 with system

    Returns:
        dict[str, np.ndarray]: A interaction matrix, r growth rate, X0 initial conditions
    """
    with h5py.File(filepath, "r") as h5:
        # 1. Read Datasets (as NumPy arrays)
        # Use [:] to load the data into memory
        data = {
            "A": h5["system/A"][:],
            "r": h5["system/r"][:],
            "X0": h5["system/X0"][:]
        }

    return data

def save_to_hdf5(
    filepath: Path,
    cfg: Config,
    deterministic: torch.Tensor,
    stochastic: torch.Tensor,
    A: torch.Tensor,
    r: torch.Tensor,
    X0: torch.Tensor,
    sigma: torch.Tensor,
    metadata: dict[str, str],
) -> None:
    """Persist simulation outputs and configuration to an HDF5 file."""
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(filepath, "w") as h5:
        h5.create_dataset("trajectories/deterministic", data=deterministic.cpu().numpy())
        h5.create_dataset("trajectories/stochastic", data=stochastic.cpu().numpy())
        h5.create_dataset("system/A", data=A.cpu().numpy())
        h5.create_dataset("system/r", data=r.cpu().numpy())
        h5.create_dataset("system/X0", data=X0.cpu().numpy())
        h5.create_dataset("system/sigma", data=sigma.cpu().numpy())

        cfg_group = h5.create_group("config")
        for key, value in vars(cfg).items():
            if value is None:
                continue
            cfg_group.attrs[key] = str(value)

        meta_group = h5.create_group("metadata")
        for key, value in metadata.items():
            meta_group.attrs[key] = value

