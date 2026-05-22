import math
import torch
import pytest

# Mark every benchmark in this file
pytestmark = pytest.mark.benchmark


from theomodels.core import simulate

# 1. We create a fixture for inputs so we aren't measuring PyTorch's allocation time
@pytest.fixture
def base_inputs():
    torch.manual_seed(42)  # Make it deterministic
    n_species = 10
    return {
        "A": torch.randn(n_species, n_species) * 0.1,
        "r": torch.randn(n_species) * 0.5,
        "X0": torch.rand(n_species) + 0.5,
        "dt": 0.01,
        "sigma": 0.05
    }

# 2. Benchmark varying 'steps' while keeping 'n_trials' constant
@pytest.mark.parametrize("steps", [10, 100, 1000])
def test_simulation_steps(benchmark, base_inputs, steps):
    
    # SETUP: We define a small inner function to pass to the benchmark fixture
    def run_simulation():
        simulate(**base_inputs, steps=steps, n_trials=1)
        
        # EXTREMELY IMPORTANT FOR PYTORCH:
        # If running on GPU, PyTorch operations are asynchronous. 
        # This line forces Python to wait until the GPU finishes calculations,
        # ensuring the benchmark timer captures the actual execution time.
        if base_inputs["X0"].is_cuda:
            torch.cuda.synchronize()

    # TEST: Benchmark runs 'run_simulation' repeatedly
    benchmark(run_simulation)


# 3. Benchmark varying 'n_trials' while keeping 'steps' constant
@pytest.mark.parametrize("n_trials", [1, 10, 100])
def test_simulation_trials(benchmark, base_inputs, n_trials):
    
    def run_simulation():
        simulate(**base_inputs, steps=100, n_trials=n_trials)
        if base_inputs["X0"].is_cuda:
            torch.cuda.synchronize()

    benchmark(run_simulation)



# 1. Map simple string labels to actual hardware checks
def get_device(device_name: str) -> torch.device:
    if device_name == "gpu":
        if not torch.cuda.is_available():
            pytest.skip("GPU requested but CUDA is not available on this machine")
        return torch.device("cuda")
    return torch.device("cpu")


# 2. Parametrize both dimensions at the test level
@pytest.mark.parametrize("n_species", [2, 4, 8, 16])
@pytest.mark.parametrize("device_type", ["cpu", "gpu"])
def test_simulation_hardware_scaling(benchmark, n_species, device_type):
    
    # SETUP: Target the correct hardware string
    device = get_device(device_type)
    
    # Generate mock tensors directly on the target hardware
    inputs = {
        "A": torch.randn(n_species, n_species, device=device) * 0.1,
        "r": torch.randn(n_species, device=device) * 0.5,
        "X0": torch.rand(n_species, device=device) + 0.5,
        "dt": 0.01,
        "steps": 100,       # Static baseline for matrix testing
        "n_trials": 5,      # Static baseline for matrix testing
        "sigma": 0.05
    }

    # EXECUTION: Define the isolated timed step
    def run_simulation():
        # Unpack the inputs dictionary into the function
        simulate(**inputs)
        
        # Don't let the GPU cheat the timer! Force sync.
        if device.type == "cuda":
            torch.cuda.synchronize()

    # TEST: Ask pytest-benchmark to evaluate it exactly 10 times per round
    benchmark.pedantic(run_simulation, iterations=1, rounds=10)