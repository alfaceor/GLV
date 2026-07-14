# GLV Project: Simulation and Inference of Generalized Lotka–Volterra Models
![Build Status](https://github.com/alfaceor/GLV/actions/workflows/ci.yml/badge.svg)
[![codecov](https://codecov.io/gh/alfaceor/GLV/branch/main/graph/badge.svg)](https://codecov.io/gh/alfaceor/GLV)
i[![Python Version](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue.svg)](https://www.python.org/)



## Generalized Lotka=Volterra (GLV)

The human gut microbiota consists of hundreds of interacting species, creating a high-dimensional system where the number of parameters in a Generalized Lotka-Volterra (GLV) model grows quadratically ($N^2$), leading to overfitting and poor identifiability.

While traditional regression methods fail to capture the non-linearities and stochastic nature of these communities, modern machine learning tools like SINDy and Neural ODEs offer a potential solution, yet their performance on ecologically constrained systems remains untested.

This project will establish a comparative framework to determine which dimensionality-reduction and inference techniques best preserve the 'biological truth' of the system while maintaining predictive accuracy under perturbation.


## Overview

This repository implements numerical simulations and analysis pipelines for **Generalized Lotka–Volterra (GLV)** models, with applications to ecological dynamics.

The project combines:

* High-performance numerical solvers (C++)
* Model and orchestration layer (Python)
* Statistical analysis and visualization (R)
* Reproducible reporting (Quarto)

The goal is to provide a **reproducible workflow** from simulation to inference analysis and  publication-quality results.

## Repository Structure


```text
├── src/                # Core Python package 
├── test/               # Unit tests
├── data/               # Input datasets (raw / external)
├── results/            # Intermediary results
├── docs/               # Project documentation and planning
├── Notebooks/          # Exploratory analysis
├── Prompts/            # AI prompts used in project
├── NumericalSolutions/ # High-performance numerical solvers (C++)
├── Synopsis/           # Quarto manuscript (methods, results, etc)
└── Bibliography/       # Bibliography used in project (only references, no pdfs)
```

---

## Installation

### 1. Clone the repository

```bash
git clone git@github.com:alfaceor/GLV.git
cd GLV
```

### 2. Set up Python environment

Create a conda environment (optional):

```bash
conda env create -f environment.yml
conda activate glv
```

Install package with pip

```bash
pip install -e .
```



### 3. Build C++ solver (optional but recommended)

<!-- FIXME: TDB if C++ implementation will be in this repo or not  -->

```bash
cd cpp
mkdir build && cd build
cmake ..
make
```

---

## Quick Start

### Theoretical simulations

To run numerical experiment

```{bash}
glv-simulate n_species=2 noise=single noise.index=0 noise.std=0.1 extft=none
```

by default, this command will generate a simulation with random  a systems with n_species=2, noise=none, extft=none


```{bash}
glv-simulate --multirun n_species=2 noise=single noise.index=0,1 noise.std=0.01,0.05,0.1,0.5 extft=none
```



```{bash}
mlflow ui
```


<!-- TODO: Future implementation of pipeline, to make reproducible using 'make' -->

Run the full pipeline:



## Reproducibility

This project is designed to be reproducible:

* All figures in the paper are generated from scripts in `scripts/`
* Intermediate results are stored in `results/`
* No manual steps are required beyond running the pipeline

If results differ across systems, please check:

* Python/R versions
* Compiler version for C++
* Random seeds (TODO: specify where controlled)


## Example Usage

Run a GLV simulation programmatically:

```bash
glv-simulate n_species=2 noise=all noise.std=0.5
```

```python
from theomodels.core import simulate

results = simulate(
    params=TODO,
    t_span=TODO
)
```

---

## Testing

Run unit tests with:

```bash
pytest -q
```

---

## Scientific Report

The full manuscript is located in:

```text
Synopsis/
```

To render:

```bash
quarto render Synopsis/
```

---

## Dependencies

* Python: numpy, scipy, pandas, matplotlib
* R: TODO (list packages)
* C++: CMake, standard library (C++11+)
* Quarto (for report generation)


## Roadmap

* [x] Add stochastic simulations
* [ ] Improve parameter inference methods
* [ ] Integrate C++ solver with Python bindings
* [x] Add CI for automated testing
* [ ] Improve documentation and examples

## Contributing

Contributions are welcome. Please:

1. Open an issue to discuss changes
2. Submit a pull request with clear description
3. See CONTRIBUTING.md file

## License

This project is licensed under the GLP v3 (see LICENSE file).

## References

See `Bibliography/GLV_references.bib` for full bibliography.

## Contact

Carlos Olivares






# Installation

```{bash}
pip install -e .
```


# Run simulations

To add noise to a single species

```{bash}
# Single species
glv-simulate noise=single noise.id=2 noise.std=0.5
```

To add noise to all species

```{bash}

# All particles
glv-simulate noise=all noise.std=0.5
```

To add stochastic noise to selected species

```{bash}

# Selected particles
glv-simulate noise=selected noise.selid='{0:0.1,3:0.9}' noise.default_std=0.5
```





## How to infer the interactions network?

Based on deterministic and stochastic simulations of a Generalized Lotka-Volterra system, determine the conditions to predict accurately the macroscopic ecological consequences of perturbations.


## How to establish uncertainty on the inferred parameters?


## How to balance interpretability and predictive power? (trade-off)


## What is the dimensionaly reduced system?

- How to encode this high-dimensional data. PCA, use autoenconders?
- simple regression
- Non-linear mixed effects
- Neural Networks
- SYNDy or MANDy?
- Graphical Neural Networks (GNN)
- Neural ODE 
- Physics Informed Neural Network (PINN)
- GPINN?

## Phase 0: Make Deterministic and Stochastic Simulations



## Phase 1: The "Methodological Comparison" Focus

**Title:** Benchmarking Sparse vs. Black-Box Inference for Gut Microbiota Stability.
**Research Objective:** To systematically compare the accuracy of **SINDy** (Sparse Identification of Non-linear Dynamics) and **Neural ODEs** in recovering the interaction coefficients ($A_{ij}$) of a gut-like microbial community from noisy, time-resolved data.
**Key Question:** Under what levels of measurement noise and sampling frequency does a physics-informed model outperform a purely data-driven neural network?


## Phase 2: The "Ecological Perturbation" Focus

**Title:** Physics-Informed Neural Networks (PINNs) for Predicting State Transitions in the Gut Microbiota.
**Research Objective:** To develop a hybrid modeling framework that uses **Autoencoders** for dimensionality reduction and **PINNs** to simulate the macroscopic response of the gut microbiota to external perturbations (e.g., antibiotics or dietary shifts).
**Key Question:** Can a reduced-order model effectively predict a "tipping point" or a loss of diversity without modeling every single rare species?


## Phase 3: The "Latent Space" Focus

**Title:** Uncovering the Latent Dynamics of Microbial Ecology through Manifold Learning.
**Research Objective:** To evaluate the effectiveness of **Graph Neural Networks (GNNs)** and **Latent ODEs** in encoding high-dimensional taxonomic data into a low-dimensional manifold that preserves the underlying Generalized Lotka-Volterra physics.
**Key Question:** Does a latent-space representation of the microbiota provide a more stable prediction of long-term community equilibrium than species-level GLV models?



# Code execution

## Run the full pipeline once
python simulate.py --multirun model=GLVdeterministic,GLVstochastic
python prepare_data.py ...
python train.py --multirun inference=nlme,neural_ode,GNN
python evaluate.py

## Render the report
cd Synopsis
quarto render report.qmd --to pdf

## Iterate on prose only — chunks use cached outputs
quarto render report.qmd --to pdf   # fast, no HDF5 reads

## New simulation result arrived — rerun affected chunks only
quarto render report.qmd --to pdf   # freeze:auto handles this



# Want to contribute?
See CONTRIBUTING.md
