from abc import ABC
from dataclasses import dataclass, field
from typing import Any
from hydra.core.config_store import ConfigStore


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
    noise: NoiseConfig = field(default_factory=NoneNoise)



def register_configs() -> None:
    """Register all configs in the ConfigStore."""
    cs = ConfigStore.instance()
    cs.store(name="base_config", node=Config)
    cs.store(group="noise", name="none", node=NoneNoise)
    cs.store(group="noise", name="single", node=SingleNoise)
    cs.store(group="noise", name="all", node=AllNoise)
    cs.store(group="noise", name="selid", node=SelIdNoise)


def register_configs_buildsystem() -> None:
    cs = ConfigStore.instance()
    cs.store(name="base_config", node=Config)
