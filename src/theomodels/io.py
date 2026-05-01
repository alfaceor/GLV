import torch
# from dataclasses import dataclass, field
# from typing import Any
from theomodels.config import Config
import h5py
import numpy as np
from pathlib import Path
from typing import Union
import platform
import numpy as np
from typing import TypedDict



class GLVParams(TypedDict):
    A: np.ndarray
    r: np.ndarray
    X0: np.ndarray

class NoiseParams(TypedDict):
    sigma: np.ndarray


_HDF5_KEY_A   = "system/A"
_HDF5_KEY_R   = "system/r"
_HDF5_KEY_X0  = "system/X0"
_HDF5_KEY_TRAJ = "traj"
_HDF5_KEY_TIME = "time"
_HDF5_KEY_SIGMA = "noise/sigma"



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
    sigma = torch.full((cfg.n_species,), cfg.noise.std, device=device)
    return sigma

def build_sigma_single(cfg: Config, device: torch.device) -> torch.Tensor:
    sigma = torch.zeros(cfg.n_species, device=device)
    sigma[cfg.index] = cfg.std
    return sigma
    # raise NotImplementedError


def build_sigma_selid(cfg: Config, device: torch.device) -> torch.Tensor:    
    base_std = cfg.noise.default_std
    sigma = torch.full((cfg.n_species,), base_std, device=device)
    # sigma = torch.full((cfg.n_species,), cfg.default_std, device=device)
    noise_map = parse_noise_map(cfg.noise_map_str, cfg.n_species)
    try:
        for key, value in noise_map.items():
            sigma[key] = value
    except Exception as e:
        # print(e)
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
    elif cfg.noise._target_ == "SingleNoise":
        sigma = build_sigma_single(cfg, device)
    elif cfg.noise._target_ == "SelIdNoise":
        sigma = build_sigma_selid(cfg, device)
    elif cfg.noise._target_ == "NoneNoise":
        sigma = None
    else:
        raise ValueError("Unknow noise configuration")
    return sigma


def save_system_to_hdf5(
    filepath: Union[Path, str],
    A: np.ndarray,
    r: np.ndarray,
    X0: np.ndarray,
) -> None:
    """ Save growth rate, interaction matrices with implicit carrying capacities
    Args:
        filepath (Path): filepath to save data in h5 file format
        A (np.ndarray): Interaction matrix between species
        r (np.ndarray): growth rate
        X0 (np.ndarray): initial conditions
    """
    if isinstance(filepath, str):
        filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    if np.any(np.isnan(A)) or np.any(np.isinf(A)):
        raise ValueError("Interaction matrix A contains NaN or Inf values")

    if A.shape[0] != A.shape[1]:
        raise ValueError(f"A must be square, got shape {A.shape}")
    if A.shape[0] != r.shape[0]:
        raise ValueError(f"A and r dimensions must match: {A.shape[0]} vs {r.shape[0]}")

    with h5py.File(filepath, "w") as h5:
        dtattr = h5py.string_dtype(encoding="utf-8")
        h5.create_dataset(_HDF5_KEY_A, data=A, dtype=np.float64, compression="gzip", compression_opts=4)
        h5.create_dataset(_HDF5_KEY_R, data=r, dtype=np.float64, compression="gzip", compression_opts=4)
        h5.create_dataset(_HDF5_KEY_X0, data=X0, dtype=np.float64, compression="gzip", compression_opts=4)
        # FIXME: Evaluate if it is better to create a group and place it as metadata all attributes bellow
        h5.attrs.create("torch_version", torch.__version__, dtype=dtattr)
        # h5.attrs["torch_version"] = np.string_(torch.__version__) #.encode("utf-8")
        h5.attrs["numpy_version"] = np.__version__
        h5.attrs["python_version"] = platform.python_version()
        h5.attrs["schema_version"] = "1.0"
        h5.attrs["created_with"] = "theomodels"


def load_system_from_hdf5(
    filepath: Union[Path, str]
) -> GLVParams:
    """Load system data from hdf5
    Args:
        filepath (Path): Path of hdf5 file
    Returns:
        dict: dictionary of numpy arrays with keys
        A: interaction matrix between species,
        r: growth rate,
        X0: initial conditions
    """
    if isinstance(filepath, str):
        filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"HDF5 file not found: {filepath}")
    with h5py.File(filepath, "r") as h5:
        data = {
            "A": h5[_HDF5_KEY_A][:],
            "r": h5[_HDF5_KEY_R][:],
            "X0": h5[_HDF5_KEY_X0][:]
        }
    return data


def save_traj_to_hdf5(
    filepath: Path,
    cfg: Config,
    traj: np.ndarray,
    time: np.ndarray,
    A: np.ndarray,
    r: np.ndarray,
    X0: np.ndarray,
    sigma: Union[np.ndarray, None],
    metadata: dict[str, str],
) -> None:
    """Save trajectory data to hdf5
    Args:
        filepath (Path): filepath to save data in h5 file format
        cfg (Config): configuration object
        traj (np.ndarray): trajectory data
        time (np.ndarray): time points
        A (np.ndarray): Interaction matrix between species
        r (np.ndarray): growth rate
        X0 (np.ndarray): initial conditions
        sigma (np.ndarray or None): noise term, or None if deterministic
        metadata (dict[str, str]): dictionary of metadata key-value pairs
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(filepath, "w") as h5:
        h5.create_dataset(_HDF5_KEY_TRAJ, data=traj)
        h5.create_dataset(_HDF5_KEY_TIME, data=time)
        h5.create_dataset(_HDF5_KEY_A, data=A)
        h5.create_dataset(_HDF5_KEY_R, data=r)
        h5.create_dataset(_HDF5_KEY_X0, data=X0)
        if sigma is not None:
            h5.create_dataset(_HDF5_KEY_SIGMA, data=sigma)
        cfg_group = h5.create_group("config")
        for key, value in vars(cfg).items():
            if value is None:
                continue
            cfg_group.attrs[key] = str(value)
        meta_group = h5.create_group("metadata")
        for key, value in metadata.items():
            meta_group.attrs[key] = value

def explore(name, obj):
    """Print hdf5 Dataset attributes, shape and dtype

    Args:
        name (_type_): _description_
        obj (_type_): _description_
    """
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

