# Package index

## Package

Package overview

- [`bhfvar-package`](https://joonho112.github.io/bhfvar/reference/bhfvar-package.md)
  [`bhfvar`](https://joonho112.github.io/bhfvar/reference/bhfvar-package.md)
  : bhfvar: Bayesian Hybrid Framework for Variance Decomposition

## Model Compilation

Functions for compiling the Stan model

- [`compile_bhf_model()`](https://joonho112.github.io/bhfvar/reference/compile_bhf_model.md)
  : Compile the BHF Stan Model
- [`get_stan_file_path()`](https://joonho112.github.io/bhfvar/reference/get_stan_file_path.md)
  : Get Path to Stan Model File

## Data Preparation

Functions for preparing survey data

- [`prepare_bhf_data()`](https://joonho112.github.io/bhfvar/reference/prepare_bhf_data.md)
  : Prepare Data for BHF Model

## Model Fitting

Functions for fitting the Bayesian model

- [`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md)
  : Fit the Bayesian Hybrid Framework Model

## Result Extraction

Functions for extracting and summarizing results

- [`variance_decomposition()`](https://joonho112.github.io/bhfvar/reference/variance_decomposition.md)
  : Extract Article-Aligned Variance Decomposition Results
- [`domain_estimates()`](https://joonho112.github.io/bhfvar/reference/domain_estimates.md)
  : Extract Domain-Level Estimand A or B Probabilities
- [`overall_estimate()`](https://joonho112.github.io/bhfvar/reference/overall_estimate.md)
  : Extract an Overall Estimand A or B Population Probability
- [`log_lik()`](https://joonho112.github.io/bhfvar/reference/log_lik.md)
  : Extract Raw or Pseudo Pointwise Log-Likelihood Contributions

## Plots

Graphics for the recommended workflow (needs ggplot2)

- [`bhf_plot_variance()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_variance.md)
  : Plot the Probability-Scale Variance Decomposition
- [`bhf_plot_icc()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_icc.md)
  : Plot Latent-Scale Intraclass Correlations
- [`bhf_plot_domains()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_domains.md)
  : Plot Domain Estimates
- [`bhf_plot_shrinkage()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_shrinkage.md)
  : Plot Shrinkage of Domain Estimates
- [`bhf_plot_astar_sensitivity()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_astar_sensitivity.md)
  : Plot the A\* Sampling-Variance Sensitivity

## Diagnostics

Sensitivity of A\* to its fixed inputs

- [`bhf_astar_sensitivity()`](https://joonho112.github.io/bhfvar/reference/bhf_astar_sensitivity.md)
  : Sensitivity of Estimand A\* to the Sampling Variances

## Utilities

Helper functions

- [`calc_eff_n()`](https://joonho112.github.io/bhfvar/reference/calc_eff_n.md)
  : Calculate Effective Sample Size

## Data

Example datasets

- [`bhf_synthetic_data`](https://joonho112.github.io/bhfvar/reference/bhf_synthetic_data.md)
  : Synthetic Survey Data for BHF Package Examples
