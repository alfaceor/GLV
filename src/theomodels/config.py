from abc import ABC
from dataclasses import dataclass, field
from typing import Any
from hydra.core.config_store import ConfigStore

from omegaconf import MISSING


# TODO: Add clases to manage the external perturbation
# Initiate with 2 cases:
# - No external perturbation.
# - Constant value perturbation during a fixed time.
# - Add sin or cos perturbation or periodic in general.
# Like the noise the perturbation should be over one species or multiple and with different amplitudes.
@dataclass
class ExtForceConfig:
    pass

@dataclass
class EFNone(ExtForceConfig):
    type: str = "EFNone"

@dataclass
class EFSingle(ExtForceConfig):
    f0: float = 0.0
    # psteps: int | None = None
    psteps: int = 1
    index : int = 0
    type: str = "EFSingle"

@dataclass
class EFAll(ExtForceConfig):
    f0: float = 0.0
    # psteps: int | None = None # perturbation steps. None should be all time steps
    psteps: int = 1
    type: str = "EFAll"

@dataclass
class EFSelId(ExtForceConfig):
    map_str: str = "0:0.0"
    defval: float = 0.0 # default value
    psteps: int = 1
    type: str = "EFSelId"


@dataclass
class NoiseConfig:
    pass  # Base class for type-hinting

@dataclass
class NoneNoise(NoiseConfig):
    type: str = "NoneNoise"

@dataclass
class SingleNoise(NoiseConfig):
    std: float = 0.1
    index: int = 0
    type: str = "SingleNoise"
    # hydra cli = single

@dataclass
class AllNoise(NoiseConfig):
    std: float = 0.1
    type: str = "AllNoise" # Helper tag for logic
    # hydra cli = all

@dataclass
class SelIdNoise(NoiseConfig):
    noise_map_str: str = "0:0.0"        # e.g. "0:0.1;2:0.4"
    # This allows specific overrides
    # noise_map: dict[int, float] = field(default_factory=dict)
    default_std: float = 0.0
    type: str = "SelIdNoise"
    # hydra cli = selected


# TODO: Implement system parameters to select different base models
@dataclass
class SystemConfig(ABC):
    pass # Abstract Base Class for system configuration

# FIXME: TO BE DECIDED
# NOTE: Future implementation
# TBD: To be decided if the class will use a file or random generator
# meanwhile defined as an abstract class
@dataclass
class GLVSystem(ABC):
    h5input: str = "<path>.h5"
    mode: str = "feasible"
    type: str = "GLVSystem"


# FIXME: TO BE DECIDED
# NOTE: Future implementation
@dataclass
class CRSystem(ABC):
    h5input: str = "<path>.h5"
    type: str = "CRSystem"


@dataclass
class Config:
    defaults: list[Any] = field(default_factory=lambda: [
        "_self_",
        {"noise": "none"},   # declares noise as a swappable group
        {"extft": "none"}
    ])
    n_species: int = 4
    dt: float = 0.01
    steps: int = 50
    n_trials: int = 10
    # noise_std: float = 0.1
    seed: int | None = 42
    system_mode: str = "feasible" # FIXME: Change system mode to system type
    device: str = "auto"
    save_results: bool = True
    output_dir: str = "data/simul"
    run_name: str = "glv_sim"
    noise: NoiseConfig = field(default_factory=NoiseConfig)
    extft: ExtForceConfig = field(default_factory=ExtForceConfig)


def register_configs() -> None:
    """Register all configs in the ConfigStore."""
    cs = ConfigStore.instance()
    cs.store(name="base_config", node=Config)
    cs.store(group="noise", name="none", node=NoneNoise)
    cs.store(group="noise", name="single", node=SingleNoise)
    cs.store(group="noise", name="all", node=AllNoise)
    cs.store(group="noise", name="selid", node=SelIdNoise)
    cs.store(group="extft", name="none", node=EFNone)
    cs.store(group="extft", name="single", node=EFSingle)
    cs.store(group="extft", name="all", node=EFAll)
    cs.store(group="extft", name="selid", node=EFSelId)


def register_configs_buildsystem() -> None:
    cs = ConfigStore.instance()
    cs.store(name="base_config", node=Config)
