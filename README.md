# Honest Adaptive Energy Testing for Equality of Multivariate Distributions: A Partitioned Inference Protocol

This repository provides the complete, reproducible R implementation and replication code for the paper:
> Salehi, M. (2026) Honest Adaptive Energy Testing for Equality of Multivariate Distributions: A Partitioned Inference Protocol.

---

## Overview of the Methodology
Distance-based energy statistics—specifically the Distance Components (DISCO) framework—depend fundamentally on the distance exponent alpha in (0, 2]. Selecting this structural parameter adaptively on the full sample breaks label exchangeability under the sharp null hypothesis, leading to severe Type I error inflation.

This repository provides:
1. F_{E,\mathrm{split}}$ Protocol: The foundational single-split permutation test guaranteeing finite-sample validity under $H_0$.
2. $F_{E,\mathrm{multi}}$ Protocol: The multi-split aggregation framework combining split-specific $p$-values via the Cauchy combination test (Liu & Xie, 2020), which eliminates partition variability and recovers over 95% of testing efficiency.

---

## Repository Structure
- `run_definitive_honest_energy_simulations.R`: Complete factorial Monte Carlo simulation study across dimensions d in {10, 50}, sample sizes  in {30, 60, 100}, and diverse distributional geometries (Gaussian, Student's t_3, and heavy-tailed Cauchy).
- `run_real_golub_multisplit.R`: Multi-split sensitivity analysis on the high-dimensional Golub leukemia microarray dataset (p = 3051, N = 38).
- `run_multisplit_case_study.R`: Multi-split sensitivity analysis on the multi-group oncology cohort (k = 3, N = 150, p = 100).
- `definitive_simulation_results.csv`: Precomputed simulation results across 24 factorial scenarios (1000 replications, 500 permutations, 50 splits per run).

---

## Dependencies
The code requires R (version > 4.2.0) and the following CRAN packages:

install.packages(c("tidyverse", "energy", "mvtnorm", "furrr", "future", "parallel"))

---

## Reproducibility
All scripts utilize the L'Ecuyer-CMRG random number generator with fixed seeds (MAIN_SEED <- 110) to guarantee deterministic cross-platform reproducibility.
```R
