# Online inference for dynamic categorical regression

This repository contains the anonymized code and numerical outputs accompanying the manuscript.

## Repository structure

- `scripts/` contains the filtering, smoothing, simulation, and EM routines.
- `experiments/` contains the scripts used for the simulation studies and the Expedia application. The `.com` files are generic SLURM submission templates.
- `data/` contains the processed Expedia data and the saved numerical outputs used to construct the tables and figures.
- `figures/` contains the figures included in the manuscript and Supplementary Materials.

## Reproducing tables and figures from saved output

Run the following commands from the repository root:

```sh
Rscript filtering_tables.R
Rscript filtering_traj.R
Rscript repeated_EM.R
Rscript repeated_EM_convergence.R
Rscript expedia_inference.R
```

The scripts reproduce, respectively:

- the filtering-accuracy and computational-cost table;
- the three filtering-trajectory figures in the Supplementary Materials;
- the covariance-recovery table and figure;
- the EM convergence diagnostic in the Supplementary Materials;
- the Expedia EM diagnostic, estimated covariance matrix, and dynamic-coefficient figure.

`expedia_inference.R` reuses the saved 500-iteration EM fit when it is available. If the saved fit is absent, the script reruns the estimation before creating the figures.

## Reproducing the numerical experiments

The scripts in `experiments/` generate the saved output consumed by the scripts above:

- `filtering_likelihood_T10_m10_parallel.R`, `filtering_likelihood_T50_m10_parallel.R`, and `filtering_likelihood_T50_m50_parallel.R` run the filtering comparison.
- `filtering_traj_likelihood_T10_m10.R`, `filtering_traj_likelihood_T50_m10.R`, and `filtering_traj_likelihood_T50_m50.R` generate one representative dataset for each filtering-trajectory figure.
- `repeated_EM_Sigma_parallel.R` runs one covariance-recovery regime; pass `10`, `30`, or `60` as its command-line argument.
- `expedia_inference.com` is a SLURM template that runs the 500-iteration Expedia EM fit on the first 50 days through the root-level `expedia_inference.R` script.

The filtering comparison uses 96 independent state and response replications with a common covariate design and realized sample-size sequence. The covariance-recovery study uses 48 replications. Parallelization is across replications only.

## Software requirements

The code is written in R and uses the following contributed packages:

`MASS`, `abind`, `mvtnorm`, `TruncatedNormal`, `tlrmvnmvt`, `magic`, `matrixcalc`, `future`, `future.apply`, `dplyr`, `tidyr`, `transport`, `expm`, `ggplot2`, `patchwork`, and `xtable`.

All random-number seeds and simulation settings are specified in the scripts.
