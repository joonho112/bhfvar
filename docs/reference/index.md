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
- [`validate_input_data()`](https://joonho112.github.io/bhfvar/reference/validate_input_data.md)
  : Validate Input Data
- [`validate_stan_data()`](https://joonho112.github.io/bhfvar/reference/validate_stan_data.md)
  : Validate Stan Data

## Model Fitting

Functions for fitting the Bayesian model

- [`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md)
  : Fit the Bayesian Hybrid Framework Model

## Result Extraction

Functions for extracting and summarizing results

- [`variance_decomposition()`](https://joonho112.github.io/bhfvar/reference/variance_decomposition.md)
  : Extract Variance Decomposition Results
- [`domain_estimates()`](https://joonho112.github.io/bhfvar/reference/domain_estimates.md)
  : Extract Domain-Level Estimates
- [`overall_estimate()`](https://joonho112.github.io/bhfvar/reference/overall_estimate.md)
  : Extract Overall Population Estimate
- [`log_lik()`](https://joonho112.github.io/bhfvar/reference/log_lik.md)
  : Extract Log-Likelihood for LOO-CV

## Utilities

Helper functions

- [`calc_eff_n()`](https://joonho112.github.io/bhfvar/reference/calc_eff_n.md)
  : Calculate Effective Sample Size

## Data

Example datasets

- [`bhf_synthetic_data`](https://joonho112.github.io/bhfvar/reference/bhf_synthetic_data.md)
  : Synthetic Survey Data for BHF Package Examples
