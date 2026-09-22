# Empirical study

This directory contains the data, implementation, estimation scripts, and
saved results for the empirical analysis in **Shrinkage estimation and
inference for network-linked data: a robust spectral embedding approach**.
The empirical network is constructed from the 25-nearest-neighbor information
stored in `data/DataList.rda`.

## Directory contents

- `data/DataList.rda`: processed empirical data. The object `dataList`
  contains the response and covariates in `X.data` and the network information
  in `Wmatrices$knn25`.
- `functions/`: implementations of BAR-SP and the comparison methods, together
  with spectral-projection and utility functions.
- `APP_stat.R`: main empirical estimation with the network dimension fixed at
  `K = 25`. It also runs the network-perturbation experiment.
- `APP_stat ic.R`: alternative empirical estimation that selects `K` by an
  information criterion over `K = 0, ..., 50` and saves the complete criterion
  path and diagnostic plot.
- `APP_analysis.R`: post-processing and figure-generation code based on the
  saved empirical estimates. Its default input is `results/hete`.
- `results/homo/` and `results/hete/`: saved results from the fixed-`K`
  analysis, separated by the homoskedasticity diagnostic.
- `IC_K_selection/<homo|hete>/n = <n>, K = <K>/`: saved results from the
  information-criterion analysis. The current run has `n = 255`, is classified
  as heteroskedastic, and selects `K = 23`.
- `disturb agg.csv`: aggregated network-perturbation results used by the
  empirical analysis.

## R dependencies

The estimation scripts use the following packages:

```r
install.packages(c(
  "RSpectra", "readr", "nleqslv", "optimx", "Matrix", "ggplot2"
))
```

`APP_analysis.R` additionally uses `dplyr`, `maps`, `ggrepel`, `viridis`, and
`ggforce`.

## Reproducing the analysis

Run the scripts with this directory as the R working directory, because they
load `functions/` and `data/` using relative paths.

```r
setwd("path/to/BAR-SP/empirical study")

# Fixed K = 25 analysis and network-perturbation experiment
source("APP_stat.R")

# Information-criterion selection of K
source("APP_stat ic.R")

# Figures and post-processing for the fixed-K results
source("APP_analysis.R")
```

Both estimation scripts set the random seed to `123`. Before running
`APP_analysis.R` on a different result set, update its `results_prefix` value
to the desired output directory.

## Main outputs

The output location is determined by the heteroskedasticity diagnostic. The
information-criterion script further records the sample size and selected `K`
in the directory name.

- `metric.csv`: summary metrics for BAR-SP, SP, OLS, and SIM. Its six rows are
  the network-effect p-value, estimated rank, active-coefficient ratio,
  selected regularization parameter, out-of-sample MSE, and in-sample MSE.
- `param10.csv`: selected coefficient estimates, p-values, variable indices,
  and feature names for the four methods.
- `proj measure.csv`: projection measures returned by BAR-SP.
- `alpha.csv`, `xi.csv`, and `BAR residual.csv`: estimated network effects,
  latent effects, and BAR-SP residuals.
- `BAR not sig SP sig.csv` and `BAR sig SP not sig.csv`: variables whose
  significance conclusions differ between BAR-SP and SP.
- `group effect.csv` and the `high/low network effect*.csv` files: summaries
  used to interpret coefficient grouping and network effects.
- `disturb.csv`, `disturb SP.csv`, and `disturb SIM.csv`: robustness results
  from network perturbations produced by `APP_stat.R`.
- `IC_K_selection.csv` and `IC_K_selection.pdf`: the information-criterion
  values and selection plot produced by `APP_stat ic.R`.

The saved outputs are included for reproducibility; rerunning the scripts may
take time because estimation and cross-validation are performed repeatedly.
