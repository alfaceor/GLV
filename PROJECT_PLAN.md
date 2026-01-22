# Project Plan: Generalized Lotka-Volterra to model Gut Microbiome perturbations and nutrition

## Objectives

* What are the common microbiome composition and distribution of the gut microbiome.
* One-health compare distributions with gut, skin, mouth, environment.
* Explore the connection between nutrition and the gut microbiome using computational and theoretical approaches.
* Make a review about the state of the art.
* Develop and analyze mathematical and statistical models to capture microbial dynamics and nutrient interactions.
* Use the project as a structured learning process in systems biology, statistical physics, and computational modeling.

### Non-Goals

* Not aiming to produce clinical/medical conclusions YET (focus is theoretical and computational).
* Not focused on wet-lab experimentation YET. Save it for later.

---

## Core Requirements

* **Programming Language:** Python (preferred for scientific computing, modeling, and simulations).
* **Libraries/Tools:**

  * `numpy`, `JAX`, `pytorch`, `scipy` -> numerical methods, ODE solvers, symbolic checks
  * `pandas`, `tidyverse` -> data handling
  * `matplotlib`, `seaborn`, `ggplot` -> visualization
  * `scikit-learn`, `statsmodels`, `stan` -> clustering and statistical analysis
  * `networkx` -> network representations of species interactions
  * `Rmarkdown`, `Jupyter` -> Interactive summaries or reports

---

## Architecture / Structure

Planned structure of the repository:

```
FIXME: REDEFINE IT!!!
Change the current structure of the project to get a better organization

```


---

## Methodologies

* **Clustering:** Apply clustering techniques (e.g., k-means, hierarchical, DBSCAN) to microbial abundance data.
* **Species Abundance Distributions:** Fit and analyze empirical/theoretical distributions; compare with neutral and niche models.
* **Dynamical Modeling:**

  * Generalized Lotka–Volterra (gLV) equations for microbial competition/cooperation.
  * MacArthur Consumer Resource Model for nutrient–microbe interactions.
* **Simulation & Inference:**

  * Solve ODE systems numerically.
  * Explore parameter inference using statistical physics methods.

---

## Constraints

* Code should be modular and well-documented.
* Avoid unnecessary complexity; start simple and extend.
* Ensure reproducibility (random seeds, version-controlled datasets). DVC

---

## Next Steps

1. **Setup Environment:** Create project structure, install required libraries.
2. **Lotka–Volterra Prototype:** Implement and simulate simple gLV dynamics with a small number of species.
3. **Bibliography summaries:** Organize bibliography notes to establish the key points and the state of the art.
4. **Clustering Pipeline:** Load a synthetic or public microbiome dataset and apply clustering methods.
5. **Species Abundance Analysis:** Fit and visualize distributions.
6. **Extend Models:** Add MacArthur consumer resource model simulations.
7. **Integrate Nutrition:** Introduce nutrient-resource parameters into models.

---

## Future Directions (Optional)

* Extend to stochastic differential equations (SDEs) for noise modeling.
* Explore Bayesian inference methods for parameter estimation.
* Investigate network approaches to microbiome interactions.
