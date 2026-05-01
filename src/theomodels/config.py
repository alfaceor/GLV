# config.py
from dataclasses import dataclass, field
from typing import Any
from hydra.core.config_store import ConfigStore


@dataclass
class NoiseConfig:
    pass  # Base class for type-hinting

@dataclass
class NoneNoise(NoiseConfig):
    _target_: str = "NoneNoise"

@dataclass
class SingleNoise(NoiseConfig):
    std: float = 0.1
    index: int = 0
    _target_: str = "SingleNoise"
    # hydra cli = single

@dataclass
class AllNoise(NoiseConfig):
    std: float = 0.1
    _target_: str = "AllNoise" # Helper tag for logic
    # hydra cli = all

@dataclass
class SelIdNoise(NoiseConfig):
    noise_map_str: str = "0:0.0"        # e.g. "0:0.1;2:0.4"
    # This allows specific overrides
    # noise_map: dict[int, float] = field(default_factory=dict)
    default_std: float = 0.0
    _target_: str = "SelIdNoise"
    # hydra cli = selected


@dataclass
class SystemConfig:
    pass # Abstract Base Class for system configuration



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
    system_mode: str = "feasible"
    device: str = "auto"
    save_results: bool = True
    output_dir: str = "data/simul"
    run_name: str = "glv_sim"
    noise: NoiseConfig = field(default_factory=NoiseConfig)



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
