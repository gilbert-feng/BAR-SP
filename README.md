# BAR-SP

This repository contains the R code, processed data, and saved empirical
results for **Shrinkage estimation and inference for network-linked data: a
robust spectral embedding approach**. It provides implementations of BAR-SP,
the spectral-projection estimator, and the comparison methods used in the
simulation and empirical studies.

## Repository structure

- [`simulation/`](simulation/): synthetic-data experiments for estimation,
  inference, variable selection, grouped covariates, prediction, changing
  covariate dimension, and changing sparsity. See the
  [simulation guide](simulation/README.md) for the script-to-experiment map.
- [`empirical study/`](empirical%20study/): processed empirical data,
  fixed-dimension and information-criterion analyses, saved results, and
  figure-generation code. See the
  [empirical-study guide](empirical%20study/README.md) for full details.

Each study directory contains its own `functions/` folder. Run a script with
its study directory as the R working directory so that these functions and
other relative paths resolve correctly.

## Methods included

- **BAR-SP**: the proposed shrinkage estimator with robust spectral
  projection.
- **SP**: the unpenalized spectral-projection method.
- **OLS**, **RNC**, and **SIM**: comparison methods used where applicable.

## Requirements

The code is written in R. Package requirements differ slightly across the two
studies; the complete installation commands are listed in their respective
README files. The main dependencies include `RSpectra`, `Matrix`, `optimx`,
`nleqslv`, `pracma`, `readr`, and `ggplot2`.

## Quick start

Clone the repository and run the desired script from its containing study
directory. For example:

```r
# Simulation study
setwd("path/to/BAR-SP/simulation")
source("sim_total.R")

# Empirical study
setwd("path/to/BAR-SP/empirical study")
source("APP_stat.R")
```

The scripts set random seeds internally and create their output directories
when needed. The default simulation settings use many Monte Carlo repetitions
and sample sizes as large as 4,000, so a full run can be computationally
expensive. For a quick code check, reduce `rep_time` and `n_seq` in the chosen
simulation script.

## Empirical data

The processed empirical data in `empirical study/data/DataList.rda` are based
on the dataset used by Crespo Cuaresma and Feldkircher (2013). The empirical
scripts construct a symmetric 25-nearest-neighbor network from the stored
network information.

## Reference

Crespo Cuaresma, J. and Feldkircher, M. (2013). “Spatial filtering, model
uncertainty and the speed of income convergence in Europe.” *Journal of
Applied Econometrics*, 28(4), 720–741.
