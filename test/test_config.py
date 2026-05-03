import pytest
from dataclasses import dataclass, field
from hydra import initialize, compose
from hydra.core.config_store import ConfigStore
from hydra.core.global_hydra import GlobalHydra
from omegaconf import OmegaConf
from omegaconf import DictConfig

from theomodels.config import (
    Config,
    NoiseConfig,
    NoneNoise,
    SingleNoise,
    AllNoise,
    SelIdNoise
)
from theomodels.constants import *

import json
import h5py

# ── Fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture(autouse=True)
def clear_hydra():
    GlobalHydra.instance().clear()
    yield
    GlobalHydra.instance().clear()


@pytest.fixture
def register_configs():
    cs = ConfigStore.instance()
    cs.store(name="base_config", node=Config)
    cs.store(group="noise", name="none", node=NoneNoise)
    cs.store(group="noise", name="single", node=SingleNoise)
    cs.store(group="noise", name="all", node=AllNoise)
    cs.store(group="noise", name="selid", node=SelIdNoise)


@pytest.fixture
def base_cfg(register_configs):
    with initialize(config_path=None, version_base=None):
        cfg = compose(config_name="base_config")
    return cfg


class TestHydraInstance:
    # hydra compose returns a DictConfig instance not a Config
    def test_hydra_DictConfig(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config")
        assert isinstance(cfg, DictConfig)
        assert not isinstance(cfg, Config)

    def test_hydra_conversion_to_object(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config")
        cfg_Config = OmegaConf.to_object(cfg)
        assert isinstance(cfg_Config, Config)

    def test_hydra_conversion_to_yaml(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config")
        cfg_yaml = OmegaConf.to_yaml(cfg)
        assert isinstance(cfg_yaml, str)

    def test_hydra_conversion_to_json(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config")
        cfg_json = json.dumps(
                OmegaConf.to_container(cfg, resolve=True)
            )
        assert isinstance(cfg_json, str)

        

# ── Overrides ─────────────────────────────────────────────────────────────────

class TestOverrides:

    def test_override_n_species(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["n_species=10"])
        assert cfg.n_species == 10

    def test_override_seed_to_none(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["seed=null"])
        assert cfg.seed is None

    def test_override_device(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["device=cpu"])
        assert cfg.device == "cpu"

    def test_override_noise_to_all(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["noise=all"])
        assert cfg.noise.type == "AllNoise"

    def test_override_noise_to_all_std_01(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["noise=all", "noise.std=0.5"])
        assert pytest.approx(cfg.noise.std) == 0.5

    def test_override_noise_to_single(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["noise=single"])
        assert cfg.noise.type == "SingleNoise"

    def test_override_noise_to_single_std_01(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["noise=single", "noise.std=0.5"])
        assert pytest.approx(cfg.noise.std) == 0.5

    def test_override_noise_to_single_index_01(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["noise=single", "noise.index=2"])
        assert cfg.noise.index == 2

    def test_override_noise_to_selid(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["noise=selid"])
        assert cfg.noise.type == "SelIdNoise"

    def test_override_system_mode(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["system_mode=random"])
        assert cfg.system_mode == "random"


# ── SingleNoise config group ──────────────────────────────────────

class TestSingleNoise:

    def test_noise_single_target(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["noise=single"])
        assert cfg.noise.type == "SingleNoise"

    def test_noise_single_default_std(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["noise=single"])
        assert cfg.noise.std == 0.1

    def test_noise_single_override_default_std(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(
                config_name="base_config",
                overrides=["noise=single", "noise.std=0.3"]
            )
        assert pytest.approx(cfg.noise.std) == 0.3



# ── AllNoise config group ──────────────────────────────────────

class TestAllNoise:

    def test_noise_all_target(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["noise=all"])
        assert cfg.noise.type == "AllNoise"

    def test_noise_all_default_std(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["noise=all"])
        assert cfg.noise.std == 0.1

    def test_noise_all_override_default_std(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(
                config_name="base_config",
                overrides=["noise=all", "noise.std=0.3"]
            )
        assert pytest.approx(cfg.noise.std) == 0.3


# ── SelIdNoise config group ──────────────────────────────────────

class TestSelIdNoise:

    def test_noise_selid_target(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["noise=selid"])
        assert cfg.noise.type == "SelIdNoise"

    def test_noise_selid_default_std(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(config_name="base_config", overrides=["noise=selid"])
        assert pytest.approx(cfg.noise.default_std) == 0.0

    def test_noise_selid_override_default_std(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(
                config_name="base_config",
                overrides=["noise=selid", "noise.default_std=0.3"]
            )
        assert pytest.approx(cfg.noise.default_std) == 0.3

    def test_noise_selid_default_noise_map_str(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(
                config_name="base_config",
                overrides=[
                    "noise=selid"
                ]
            )
        assert cfg.noise.noise_map_str == "0:0.0"

    def test_noise_selid_override_noise_map_str(self, register_configs):
        with initialize(config_path=None, version_base=None):
            cfg = compose(
                config_name="base_config",
                overrides=[
                    "noise=selid",
                    "noise.noise_map_str='0:0.5'"
                ]
            )
        assert cfg.noise.noise_map_str == "0:0.5"


# ── Fixtures ────────────────────────────────────────────────────────────────

@pytest.fixture
def h5_path(tmp_path):
    return tmp_path / "test_cfg.h5"


@pytest.fixture
def saved_h5(base_cfg, h5_path):
    with h5py.File(h5_path, "w") as h5:
        grp = h5.create_group("config")
        grp.attrs["yaml"] = OmegaConf.to_yaml(base_cfg)
    return h5_path


# ── Tests ───────────────────────────────────────────────────────────────────

class TestCfgHdf5:

    def test_cfg_saved(self, saved_h5):
        """Config file is created and contains expected group/attr."""
        assert saved_h5.exists()

        with h5py.File(saved_h5, "r") as h5:
            assert "config" in h5
            assert "yaml" in h5["config"].attrs

    def test_cfg_roundtrip_dictconfig(self, base_cfg, saved_h5):
        """Saved YAML can be loaded back into a DictConfig."""
        with h5py.File(saved_h5, "r") as h5:
            yaml_str = h5["config"].attrs["yaml"]

        loaded_cfg = OmegaConf.create(yaml_str)

        assert isinstance(loaded_cfg, DictConfig)

        # Stronger check: content equality
        assert OmegaConf.to_container(loaded_cfg, resolve=True) == \
               OmegaConf.to_container(base_cfg, resolve=True)

    def test_cfg_roundtrip_dataclass(self, base_cfg, saved_h5):
        with h5py.File(saved_h5, "r") as h5:
            yaml_str = h5["config"].attrs["yaml"]

        loaded_cfg = OmegaConf.create(yaml_str)
        default_cfg = OmegaConf.structured(Config)
        assert isinstance(default_cfg, DictConfig)
        # TODO: Solve issue with the the OmegaConf.merge
        # structured_cfg = OmegaConf.merge(
        #     default_cfg,
        #     loaded_cfg
        # )
        # assert isinstance(structured_cfg, DictConfig)

        # obj = OmegaConf.to_object(structured_cfg)

        # assert isinstance(obj, Config)
        # assert obj.n_species == base_cfg.n_species
