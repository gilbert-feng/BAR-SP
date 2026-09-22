# Simulation study

This directory contains the synthetic-data experiments for **Shrinkage
estimation and inference for network-linked data: a robust spectral embedding
approach**. The scripts generate stochastic block-model networks and compare
BAR-SP with the spectral-projection and other benchmark methods under
homoskedastic and heteroskedastic errors.

## Directory contents

- `functions/`: network generation, spectral projection, BAR-SP, SP, OLS,
  RNC, SIM, and evaluation utilities.
- `sim_total.R`: main Monte Carlo experiment for the one-group covariate
  setting.
- `sim_total_G.R`: main Monte Carlo experiment for the four-group covariate
  setting.
- `sim_mse_*.R`: prediction and sensitivity experiments for the one-group
  setting.
- `sim_mse_G_*.R`: corresponding experiments for the four-group setting.
- `simulation.Rproj`: optional RStudio project file.

## Experiment scripts

| Script | Setting and primary output |
| --- | --- |
| `sim_total.R` | Main one-group study: estimation, inference, support recovery, and projection measures. |
| `sim_total_G.R` | Main four-group study, additionally evaluating recovery of the coefficient-group structure. |
| `sim_mse_in sample.R` | In-sample MSE and runtime comparison across BAR-SP with sampled or fitted networks, SP, OLS, RNC, and SIM. |
| `sim_mse_out sample.R` | Out-of-sample MSE and runtime comparison for BAR-SP with sampled or fitted networks and SP. |
| `sim_mse_p.R` | BAR-SP and SP prediction error as the covariate dimension changes. |
| `sim_mse_spars.R` | BAR-SP and SP prediction error as model sparsity changes. |
| `sim_mse_G_in sample.R` | Four-group version of the in-sample prediction experiment. |
| `sim_mse_G_out sample.R` | Four-group version of the out-of-sample prediction experiment. |
| `sim_mse_G_p.R` | Four-group covariate-dimension experiment. |
| `sim_mse_G_spars.R` | Four-group sparsity experiment. |

The main scripts vary sample size, network degree, observed-network treatment
(`sample` or block-model `MLE`), network-effect presence, and error variance.
Parameters for each experiment are defined near the top of its script and can
be changed before running it.

## R dependencies

Install the packages used across the simulation scripts:

```r
install.packages(c("RSpectra", "Matrix", "optimx", "nleqslv", "pracma"))
```

The main `sim_total*.R`, dimension, and sparsity scripts only attach
`RSpectra`; the additional packages are required by the in-sample and
out-of-sample comparison scripts.

## Reproducing the simulations

Run a script with this directory as the R working directory, because every
script loads the local `functions/` directory using a relative path.

```r
setwd("path/to/BAR-SP/simulation")

# Main one-group experiment
source("sim_total.R")

# Main four-group experiment
source("sim_total_G.R")

# Example prediction experiment; quote names that contain spaces
source("sim_mse_in sample.R")
```

The scripts should normally be run one at a time. Most use random seed `1234`;
`sim_mse_G_out sample.R` uses seed `123`. The main experiments default to 800
Monte Carlo replications for each configuration, while the MSE experiments
default to 100. To perform a quick smoke test, temporarily reduce `rep_time`
and `n_seq` near the top of the selected script.

## Output structure

Output directories are created automatically and encode the experiment,
network-effect setting, network treatment, and error type. Directory names
begin with `results` for one-group experiments and `results G` for four-group
experiments; `homo/` and `hete/` distinguish the two error settings.

The main experiments produce:

- `measures.csv`: projection measures and network-effect test summaries.
- `zeros.csv`: variable-selection precision and recall; the four-group version
  also includes the group-recovery measure.
- `model_1_sd_ratio.csv` and `model_2_sd_ratio.csv`: estimated-to-empirical
  standard-deviation ratios for BAR-SP and SP.
- Per-sample-size `SP_ave_beta degree_<d>.csv` files: coefficient summaries for
  each network degree.
- `true_projection.csv` in the one-group experiment: true projection measures.

The prediction scripts write `time.csv` and `model mse.csv`. The dimension and
sparsity scripts write separate MSE tables for sampled and fitted networks.
Generated result directories are not required in order to inspect or modify
the implementation.
