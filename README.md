# Biological Dynamics Simulations

Simulate four biological systems with ordinary differential equations in R: drug absorption and elimination, ligand–receptor binding, SIRS and SEIRS population dynamics.

## Quick start

Requires R and deSolve. Tested with R 4.5.2.

Run from the repository root:

```sh
Rscript install_packages.R
Rscript run.R
```

## Models

All equations, parameter values, initial conditions and time grids are in [`analysis.R`](analysis.R). No external data is required. Change the values before each `ode()` call to explore another scenario.

| Model | Description | Time unit |
|---|---|---|
| Drug absorption | Transfer from a gut compartment to the bloodstream, with elimination | Hours |
| Receptor binding | Reversible binding at a constant ligand concentration | Seconds |
| SIRS | Infection, recovery and loss of immunity | Days |
| SEIRS | Adds an exposed compartment before infection | Days |

## Output

`results/plots.pdf` contains the model trajectories. The corresponding time series are saved as `drug_absorption.csv`, `receptor_binding.csv`, `sirs.csv` and `seirs.csv` in `results/`.

With the included parameters, the drug amount in the bloodstream peaks at **2.5 hours** on the simulation grid. These are illustrative models with specified assumptions, not fitted clinical predictions.

## Acknowledgments

University of Skövde for the learning resources.
