# bhfvar: Bayesian Hybrid Framework for Variance Decomposition

The bhfvar package implements the Bayesian Hybrid Framework for variance
decomposition in complex surveys with post-hoc domains. It provides
tools for separating substantive geographic variation from design
artifacts and sampling noise.

## Key Features

- Bayesian Pseudo-Likelihood estimation for design consistency

- Hybrid generalized linear mixed models with domain and PSU effects

- Dual Estimand Framework: Policy (A/A\*) and Descriptive (B) estimands

- De-attenuation for finite-sample variance inflation correction

- Comprehensive diagnostic and visualization tools

## Main Functions

- [`compile_bhf_model`](https://joonho112.github.io/bhfvar/reference/compile_bhf_model.md):
  Compile Stan model (once per session)

- [`prepare_bhf_data`](https://joonho112.github.io/bhfvar/reference/prepare_bhf_data.md):
  Prepare data for Stan

- [`bhf_fit`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md):
  Fit the BHF model

- [`variance_decomposition`](https://joonho112.github.io/bhfvar/reference/variance_decomposition.md):
  Extract variance components

- [`domain_estimates`](https://joonho112.github.io/bhfvar/reference/domain_estimates.md):
  Extract domain-specific estimates

## Workflow

The recommended workflow is:

1.  Compile the Stan model once per R session using
    [`compile_bhf_model()`](https://joonho112.github.io/bhfvar/reference/compile_bhf_model.md)

2.  Prepare your data using
    [`prepare_bhf_data()`](https://joonho112.github.io/bhfvar/reference/prepare_bhf_data.md)

3.  Fit the model using
    [`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md)

4.  Extract results using
    [`variance_decomposition()`](https://joonho112.github.io/bhfvar/reference/variance_decomposition.md)
    and
    [`domain_estimates()`](https://joonho112.github.io/bhfvar/reference/domain_estimates.md)

## Design Philosophy

This package uses a "defensive" programming approach where the Stan
model is compiled explicitly by the user once per session, rather than
being pre-compiled during package installation. This approach:

- Avoids rstantools caching issues

- Provides clearer error messages when compilation fails

- Ensures compatibility across different R/Stan versions

- Gives users more control over the compilation process

## References

Lee, J., & Hooper, A. (2026). Disentangling signal from noise: A
Bayesian hybrid framework for variance decomposition in complex surveys
with post-hoc domains. *Mathematics*, 14(3), 512.
[doi:10.3390/math14030512](https://doi.org/10.3390/math14030512)

## See also

Useful links:

- <https://joonho112.github.io/bhfvar>

- <https://github.com/joonho112/bhfvar>

- Report bugs at <https://github.com/joonho112/bhfvar/issues>

## Author

JoonHo Lee <jlee296@ua.edu>
