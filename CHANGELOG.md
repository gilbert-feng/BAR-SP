# Changelog

## 2026-09-22

### BAR-SP boundary condition and empirical documentation

- Updated the homogeneous and heterogeneous inference branches in both
  `empirical study/functions/BAR.R` and `simulation/functions/BAR.R` so that
  the adjusted network-effect statistic is computed when `K - r >= 1`,
  including the one-dimensional boundary case.
- Expanded `empirical study/README.md` with the data and script roles, R
  dependencies, reproduction workflow, result-directory conventions, and a
  guide to the generated output files.

### Empirical study: information-criterion selection of network dimension

- Added `empirical study/APP_stat ic.R`, an empirical-analysis variant that
  selects the network dimension `K` by minimizing an information criterion
  instead of fixing `K = 25`.
- Added export of the full information-criterion path and its diagnostic plot.
- Stored the new run separately under `empirical study/IC_K_selection/`, with
  the heteroskedasticity classification, sample size, and selected dimension in
  the directory structure so that the original results remain unchanged.
- Added the output for the current data (`n = 255`, heteroskedastic case), for
  which the criterion selects `K = 23`, together with the BAR-SP, SP, OLS, and
  SIM estimates and comparison tables.
