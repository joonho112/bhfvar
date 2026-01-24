# bhfvar 0.3.0

## Initial Release

This is the initial release of the bhfvar package, implementing the Bayesian 
Hybrid Framework for variance decomposition in complex surveys with post-hoc 
domains.

### Key Features

* **Bayesian Pseudo-Likelihood**: Design-consistent inference using survey 
  weights in a Bayesian framework

* **Hybrid GLMM**: Simultaneous modeling of substantive domain effects (states) 
  and nuisance design effects (PSUs)

* **Dual Estimand Framework**: 
  - Estimand A (Policy): Logit-scale variance decomposition
  - Estimand B (Descriptive): Probability-scale variance decomposition  
  - Estimand A* (Policy Adjusted): De-attenuated variance decomposition

* **De-attenuation**: Both implicit (Bayesian shrinkage) and explicit 
  (method-of-moments) correction for finite-sample variance inflation

* **Comprehensive Output**:
  - `variance_decomposition()`: ICC across all three estimands
  - `domain_estimates()`: Domain-specific probabilities with reliability
  - `overall_estimate()`: Population-weighted overall estimate
  - `log_lik()`: Log-likelihood for LOO-CV

### Main Functions

* `compile_bhf_model()`: Compile the Stan model
* `prepare_bhf_data()`: Prepare survey data for Stan
* `bhf_fit()`: Fit the Bayesian Hybrid Framework model
* `variance_decomposition()`: Extract variance decomposition results
* `domain_estimates()`: Extract domain-specific estimates
* `overall_estimate()`: Extract overall population estimate

### Documentation

* Comprehensive vignettes covering:
  - Introduction to the package
  - Quick start guide
  - Methodology details
  - Dual estimand framework
  - Diagnostics and model checking
  - API reference

* pkgdown website configuration for online documentation

### Dependencies

* R >= 4.0.0
* rstan >= 2.26.0
* posterior >= 1.0.0
