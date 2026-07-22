import pytest
import torch
import numpy as np
from dataclasses import dataclass, field
from theomodels.io import (
    parse_noise_map,
    build_sigma_noise,
    save_system_to_hdf5,
    load_system_from_hdf5,
    save_traj_to_hdf5,
    load_traj_from_hdf5,
    build_extft_all,
    build_extft_single,
    build_extft_selid,
    build_extft
)
from theomodels.config import (
    Config, 
    NoiseConfig, 
    NoneNoise,
    SingleNoise,
    AllNoise,
    SelIdNoise,
    ExtForceConfig,
    EFNone,
    EFSingle,
    EFAll,
    EFSelId
    )

from theomodels.io import resolve_device


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def device():
    return torch.device("cpu")

# ── parse_noise_map ───────────────────────────────────────────────────────────

class TestParseNoiseMap:

    def test_basic_parsing(self):
        result = parse_noise_map("0:0.1;2:0.4", n_species=4)
        assert result == {0: 0.1, 2: 0.4}

    def test_empty_string(self):
        result = parse_noise_map("", n_species=4)
        assert result == {}

    def test_whitespace_string(self):
        result = parse_noise_map("   ", n_species=4)
        assert result == {}

    def test_single_pair(self):
        result = parse_noise_map("1:0.5", n_species=4)
        assert result == {1: 0.5}

    def test_all_species(self):
        result = parse_noise_map("0:0.1;1:0.2;2:0.3;3:0.4", n_species=4)
        assert result == {0: 0.1, 1: 0.2, 2: 0.3, 3: 0.4}

    def test_float_values(self):
        result = parse_noise_map("0:0.123456", n_species=4)
        assert pytest.approx(result[0]) == 0.123456

    def test_zero_value(self):
        result = parse_noise_map("0:0.0", n_species=4)
        assert result == {0: 0.0}

    def test_spaces_around_separator(self):
        result = parse_noise_map("0: 0.1; 2: 0.4", n_species=4)
        assert result == {0: 0.1, 2: 0.4}

    # ── Validation errors ─────────────────────────────────────────────────────

    def test_key_out_of_range_high(self):
        with pytest.raises(ValueError, match="out of range"):
            parse_noise_map("4:0.1", n_species=4)

    def test_key_negative(self):
        with pytest.raises(ValueError, match="out of range"):
            parse_noise_map("-1:0.1", n_species=4)

    def test_invalid_format_missing_colon(self):
        with pytest.raises(ValueError, match="Invalid pair"):
            parse_noise_map("0-0.1", n_species=4)

    def test_invalid_format_empty_pair(self):
        with pytest.raises(ValueError):
            parse_noise_map("0:0.1,,2:0.4", n_species=4)

    def test_invalid_value_not_float(self):
        with pytest.raises(ValueError):
            parse_noise_map("0:abc", n_species=4)

    def test_invalid_key_not_int(self):
        with pytest.raises(ValueError):
            parse_noise_map("a:0.1", n_species=4)

    def test_parse_noise_map(self):
        str_noise_map = '0:0.1;5:0.4' # into {0: 0.1, 2: 0.4}
        assert isinstance(parse_noise_map(str_noise_map, 10), dict)
        with pytest.raises(ValueError, match=r"Key .*"):
            parse_noise_map(str_noise_map, 2)
        
        with pytest.raises(ValueError, match=r"Invalid pair .*"): 
            str_noise_map = '0:0.1:0.1:0.7;5:0.4'
            parse_noise_map(str_noise_map, 10)

    

# def test_build_sigma_noise():
#     raise NotImplementedError

# def test_build_sigma_noise():
#     raise NotImplementedError

# def test_build_sigma():
#     # Config(n_species=4, noise=uniform_noise)
#     # cfg = MapNoise(default_std=0.0, noise_map_str="0:0.1,2:0.4")
#     cfg = Config(n_species=4, noise=SelectedNoise(noise_map_str="0:0.1"))
#     sigma = build_sigma(cfg, "cpu")
#     assert sigma.shape == (4,)


# def test_save_system_to_hdf5():
#     save_system_to_hdf5()


@pytest.fixture
def all_noise():
    return AllNoise(std=0.1)


@pytest.fixture
def map_noise():
    return SelIdNoise(default_std=0.0, noise_map_str="0:0.1;2:0.4")


@pytest.fixture
def base_cfg(all_noise):
    return Config(n_species=4, noise=all_noise)


# ── build_sigma ───────────────────────────────────────────────────────────────

class TestBuildSigmaUniform:

    def test_uniform_default(self, base_cfg, device):
        sigma = build_sigma_noise(base_cfg, device)
        expected = torch.full((4,), 0.1)
        assert torch.allclose(sigma, expected)

#     def test_uniform_shape(self, base_cfg, device):
#         sigma = build_sigma(base_cfg, device)
#         assert sigma.shape == (4,)

#     def test_uniform_custom_std(self, device):
#         cfg = Config(n_species=4, noise=UniformNoise(std=0.5))
#         sigma = build_sigma(cfg, device)
#         assert torch.allclose(sigma, torch.full((4,), 0.5))

#     def test_uniform_different_n_species(self, device):
#         cfg = Config(n_species=8, noise=UniformNoise(std=0.1))
#         sigma = build_sigma(cfg, device)
#         assert sigma.shape == (8,)
#         assert torch.allclose(sigma, torch.full((8,), 0.1))

#     def test_uniform_device(self, base_cfg, device):
#         sigma = build_sigma(base_cfg, device)
#         assert sigma.device.type == "cpu"


# class TestBuildSigmaMap:

#     def test_map_basic(self, device):
#         cfg = Config(n_species=4, noise=MapNoise(noise_map_str="0:0.1,2:0.4"))
#         sigma = build_sigma(cfg, device)
#         assert pytest.approx(sigma[0].item()) == 0.1
#         assert pytest.approx(sigma[1].item()) == 0.0   # default
#         assert pytest.approx(sigma[2].item()) == 0.4
#         assert pytest.approx(sigma[3].item()) == 0.0   # default

#     def test_map_default_std(self, device):
#         cfg = Config(n_species=4, noise=MapNoise(default_std=0.2, noise_map_str="0:0.5"))
#         sigma = build_sigma(cfg, device)
#         assert pytest.approx(sigma[0].item()) == 0.5   # overridden
#         assert pytest.approx(sigma[1].item()) == 0.2   # default_std
#         assert pytest.approx(sigma[2].item()) == 0.2   # default_std
#         assert pytest.approx(sigma[3].item()) == 0.2   # default_std

#     def test_map_empty_string(self, device):
#         cfg = Config(n_species=4, noise=MapNoise(default_std=0.3, noise_map_str=""))
#         sigma = build_sigma(cfg, device)
#         assert torch.allclose(sigma, torch.full((4,), 0.3))

#     def test_map_shape(self, device):
#         cfg = Config(n_species=4, noise=MapNoise(noise_map_str="0:0.1"))
#         sigma = build_sigma(cfg, device)
#         assert sigma.shape == (4,)

#     def test_map_all_species_overridden(self, device):
#         cfg = Config(n_species=4, noise=MapNoise(noise_map_str="0:0.1,1:0.2,2:0.3,3:0.4"))
#         sigma = build_sigma(cfg, device)
#         expected = torch.tensor([0.1, 0.2, 0.3, 0.4])
#         assert torch.allclose(sigma, expected)


# class TestBuildSigmaValidation:

#     def test_unknown_target_raises(self, device):
#         cfg = Config(n_species=4, noise=MapNoise(noise_map_str=""))
#         cfg.noise._target_ = "Unknown"
#         with pytest.raises(ValueError, match="Unknown noise target"):
#             build_sigma(cfg, device)

#     def test_map_key_out_of_range(self, device):
#         cfg = Config(n_species=4, noise=MapNoise(noise_map_str="9:0.1"))
#         with pytest.raises(ValueError, match="out of range"):
#             build_sigma(cfg, device)




# ── save_system_to_hdf5 ────────────────────────────────────────────────────────────

@pytest.fixture
def sample_system():
    n = 4
    return {
        "A": np.ones((n, n)),
        "r": np.ones(n),
        "X0": np.arange(n, dtype=np.float64)
    }

class TestSaveSystemHDF5:
    
    def test_sample_system(self, sample_system):
        n = 4
        A = sample_system["A"]
        r = sample_system["r"]
        X0 = sample_system["X0"]
        assert A.shape == (n, n)
        assert r.shape == (n,)
        assert X0.shape == (n,)
        assert isinstance(A, np.ndarray)

    def test_save_system_to_hdf5(self, tmp_path, sample_system):
        tmp_filepath = tmp_path / "test_hdf5.h5"
        A = sample_system["A"]
        r = sample_system["r"]
        X0 = sample_system["X0"]
        save_system_to_hdf5(tmp_filepath, A, r, X0)
        data = load_system_from_hdf5(tmp_filepath)
        assert isinstance(data, dict)
        A_loaded = data["A"]
        r_loaded = data["r"]
        X0_loaded = data["X0"]
        assert isinstance(A_loaded, np.ndarray)
        assert A_loaded.shape == A.shape
        assert r_loaded.shape == r.shape
        assert X0_loaded.shape == X0.shape
        np.testing.assert_allclose(A, A_loaded, rtol=1e-5)
        np.testing.assert_allclose(r, r_loaded, rtol=1e-5)
        np.testing.assert_allclose(X0, X0_loaded, rtol=1e-5)

    def test_load_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_system_from_hdf5(tmp_path / "does_not_exist.h5")


class TestSaveTrajHDF5:
    # TODO: Make a test to 
    # Test the corr
    def test_save_traj_to_hdf5_with_sigma_tensor(self, tmp_path):
        # --- Arrange ---
        filename = tmp_path / "test_output.h5"
        traj = np.random.normal(size=(10, 3))
        tt = np.linspace(0, 1, 10)
        A = np.eye(3)
        r = np.zeros(3)
        X0 = np.zeros(3)
        sigma = np.array([0.1, 0.2, 0.3])
        metadata = {"experiment": "test"}
        cfg = {"version": 1.0}
        
        # --- Act ---
        save_traj_to_hdf5(
            filepath=filename,
            traj=traj,
            time=tt,
            A=A,
            r=r,
            X0=X0,
            sigma=sigma,
            metadata=metadata,
            cfg=cfg
        )

        # --- Assert ---
        assert filename.exists()



class TestLoadTrajHDF5:
    def test_load_traj_raises_FileNotFounderror(self):
        with pytest.raises(FileNotFoundError):
            load_traj_from_hdf5("noexiste.h5")
    
    # TODO: Generate a test for KeyError
    # def test_load_traj_raises_keyError(self, tmp_path):
    #     with pytest.raises(KeyError):
    #         load_traj_from_hdf5(tmp_path)

    def test_load_traj(self, sample_system, tmp_path):
        # --- Arrange ---
        filename = tmp_path / "test_output_02.h5"
        n_trials = 10
        n_tsteps = 20
        n_species = 4
        traj = np.random.normal(size=(n_trials, n_tsteps, n_species))
        tt = np.linspace(0, 1, n_tsteps)
        A = np.eye(n_species)
        r = np.ones(n_species)
        X0 = np.ones(n_species)
        sigma = np.random.normal(size=(n_tsteps,))
        metadata = {"experiment": "test"}
        cfg = {"version": 1.0}

        # --- Act ---
        save_traj_to_hdf5(
            filepath=filename,
            traj=traj,
            time=tt,
            A=A,
            r=r,
            X0=X0,
            sigma=sigma,
            metadata=metadata,
            cfg=cfg
        )

        assert filename.exists()
        data = load_traj_from_hdf5(filepath=filename)
        assert isinstance(data, dict)
        assert data["traj"].shape == (n_trials, n_tsteps, n_species)
        np.testing.assert_allclose(A, data["system"]["A"], rtol=1e-5)
        np.testing.assert_allclose(r, data["system"]["r"], rtol=1e-5)
        np.testing.assert_allclose(X0, data["system"]["X0"], rtol=1e-5)
        np.testing.assert_allclose(traj, data["traj"], rtol=1e-5)
        np.testing.assert_allclose(tt, data["time"], rtol=1e-5)

# ── resolve_device ────────────────────────────────────────────────────────────

class TestResolveDevice:

    def test_cpu(self):
        device = resolve_device("cpu")
        assert device == torch.device("cpu")

    def test_auto_returns_device(self):
        device = resolve_device("auto")
        assert device.type in ("cpu", "cuda")

    @pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
    def test_cuda(self):
        device = resolve_device("cuda")
        assert device == torch.device("cuda")


###########################3


@pytest.fixture
def extft_all_default():
    return EFAll()

@pytest.fixture
def cfg_extft_all_default(extft_all_default):
    return Config(extft=extft_all_default)

@pytest.fixture
def cfg_extft_all_newvalue():
    return Config(extft=EFAll(f0=2.0))

@pytest.fixture
def extft_single_default():
    return EFSingle()

@pytest.fixture
def cfg_extft_single_default(extft_single_default):
    return Config(extft=extft_single_default)

@pytest.fixture
def extft_selid_default():
    return EFSelId()

@pytest.fixture
def cfg_extft_selid_default(extft_selid_default):
    return Config(extft=extft_selid_default)

class TestExtFT:
    def test_build_extft_all(self, cfg_extft_all_default):
        cfg = cfg_extft_all_default
        assert isinstance(cfg, Config)
        assert isinstance(cfg.extft, ExtForceConfig)
        assert isinstance(cfg.extft, EFAll)
        assert cfg.extft.type == "EFAll"
        assert cfg.extft.f0 == 1.0
        assert cfg.extft.psteps == 1
        assert cfg.n_species == 2
        assert cfg.steps == 50
        device = "cpu"
        f_t = build_extft_all(cfg, device=device)
        assert isinstance(f_t, torch.Tensor)
        assert f_t.shape == (cfg.steps, cfg.n_species)
        assert f_t[0, 0] == cfg.extft.f0
    
    def test_build_extft_all_newvalue(self, cfg_extft_all_newvalue):
        cfg = cfg_extft_all_newvalue
        device = "cpu"
        f_t = build_extft_all(cfg, device=device)
        assert cfg.extft.f0 == 2.0
        assert f_t.shape == (cfg.steps, cfg.n_species)
        tmp_f0 = cfg.extft.f0*torch.ones(cfg.n_species, device=device)
        # assert f_t[0, 0] == cfg.extft.f0
        torch.testing.assert_close(f_t[0], tmp_f0)

    def test_build_extft_single_default(self, cfg_extft_single_default):
        cfg = cfg_extft_single_default
        assert isinstance(cfg, Config)
        assert isinstance(cfg.extft, ExtForceConfig)
        assert isinstance(cfg.extft, EFSingle)
        assert cfg.extft.type == "EFSingle"
        assert cfg.extft.f0 == 1.0
        assert cfg.extft.psteps == 1
        assert cfg.n_species == 2
        assert cfg.steps == 50
        device = "cpu"
        f_t = build_extft_single(cfg, device=device)
        assert isinstance(f_t, torch.Tensor)
        assert f_t.shape == (cfg.steps, cfg.n_species)
        assert f_t[0, 0] == cfg.extft.f0

    def test_build_extft_selid_default(self, cfg_extft_selid_default):
        cfg = cfg_extft_selid_default
        assert isinstance(cfg, Config)
        assert isinstance(cfg.extft, ExtForceConfig)
        assert isinstance(cfg.extft, EFSelId)
        assert cfg.extft.type == "EFSelId"
        assert cfg.extft.map_str == "0:0.0"
        assert cfg.extft.psteps == 1
        assert cfg.extft.defval == 0.0
        device = "cpu"
        f0 = torch.full((cfg.n_species, ), cfg.extft.defval, device=device)
        fa = torch.full((cfg.n_species, ), 0.0, device=device)
        f_t = torch.zeros((cfg.steps, cfg.n_species), device=device)
        assert isinstance(f0, torch.Tensor)
        extft_map = parse_noise_map(cfg.extft.map_str, cfg.n_species)
        assert isinstance(extft_map, dict)
        for key, value in extft_map.items():
            f0[key] = value
        # Broadcast f_0 values
        f_t[0:cfg.extft.psteps, :] = f0
        torch.testing.assert_close(f_t[0], f0)
        torch.testing.assert_close(f_t[2], fa)
        f_t2 = build_extft_selid(cfg, device=device)
        torch.testing.assert_close(f_t, f_t2)

