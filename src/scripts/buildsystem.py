from theomodels.core import (
    random_system
)

import hydra
from hydra.core.config_store import ConfigStore


# TODO:
# - Generate a script to create random feasible systems from determined conditions like the number of species, stability conditions, etc
# - Save system in a hdf5 file, preferable, csv, sqlite, etc

from theomodels.config import SystemConfig

@hydra.main(version_base=None, config_name="base_config", )
def main(cfg:Config) -> None:
    raise NotImplemented