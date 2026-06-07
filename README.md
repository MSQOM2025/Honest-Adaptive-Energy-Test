# Honest-Adaptive-Energy-Test
R implementation of the Split-Group Energy Test. It corrects selection-induced bias in non-parametric k-sample testing through an honest, sample-splitting inference framework.


# Honest Adaptive Energy Testing

This repository contains the official R implementation and simulation environment for the paper:  
**"Honest Adaptive Energy Testing: A Split-Group Framework to Correct Selection Bias in Multi-Sample Comparisons"**

## Overview
Standard non-parametric energy statistics often require tuning parameters (such as the distance exponent $\alpha$) to maximize power. However, selecting these parameters on a full dataset leads to a significant "Selection Bias," inflating Type I error rates.

Our proposed **Split-Group Energy Test ($F_{E,split}$)** provides a rigorous framework for "Honest Inference" by decoupling parameter selection (on a discovery set) from significance testing (on an independent inference set).

## Key Results
- Selection Bias Correction: We demonstrate that naive adaptive testing inflates Type I error to ~11.7%, while our framework maintains the nominal 5% level.
- Power Recovery: Our method recovers >80% of the statistical power in moderate-to-large sample regimes ($n \ge 150$).

## File Structure
- `simulation_honest_adaptive_tuning.R`: Stress-test simulation for $n=150$ with adaptive $\alpha$ selection.
- `simulation_benchmark_fixed_power.R`: Benchmark simulations for $n \in \{30, 60, 100\}$ with fixed $\alpha=1$.
- `make_plots_fe_split.R`: R script to reproduce the figures presented in the paper.
- `/results`: Contains `.csv` output files from the Monte Carlo simulations.

## Installation & Requirements
The code is written in R and requires the following packages:
```R
install.packages(c("tidyverse", "energy", "furrr", "future"))
