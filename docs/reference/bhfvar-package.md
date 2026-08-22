# bhfvar: Bayesian Hybrid Framework for Variance Decomposition

The bhfvar package implements an article-informed crossed state/stratum
model with PSUs nested within strata, pseudo-likelihood survey weights,
and probability-scale variance decompositions for post-hoc domains.

## Supported workflow

- Crossed state and stratum effects with PSU effects nested in stratum

- Globally mean-one likelihood-weight scaling (legacy D2 is deprecated)

- Probability-scale Estimands A, A\*, and B with explicit gaps

- Fixed supplied or Taylor-linearized sampling variances for A\*

- Versioned prepared-data, fit, and extractor result contracts

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

## What has and has not been established

The implementation is verified against the main-text equations and the
Stan program in Appendix B of the accompanying article: a frozen
reference oracle, algebraic property tests, and centering-invariant
tests all pass within their tested contracts.

Interval calibration is a separate question. In 400 fits from a reduced
article-aligned balanced synthetic design, pooled primary 90\\ 0.900
(95\\ ranged from 0.860 to 0.930. This met a preregistered +/-5
percentage-point tolerance; it does not establish coverage for every
condition or other designs. No variance adjustment is implemented. The
restricted-data application reported in the article is not reproduced by
this package.

## Interpretation boundary

A is a design-effect-zero reference standardization, not a causal or
policy intervention. The A-versus-B gap is a standardization-sensitivity
contrast, not an informative-sampling-specific diagnostic. Latent ICCs
use the pre-centering model scale parameters rather than empirical
centered-effect variances.

A\* conditions on supplied or Taylor-estimated sampling variances; their
estimation uncertainty is not propagated. Pseudo-posterior intervals and
pseudo log likelihood do not automatically have ordinary posterior
coverage or ordinary observation-level LOO interpretations.

## Compilation design

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

- <https://joonho112.github.io/bhfvar/>

- <https://github.com/joonho112/bhfvar>

- Report bugs at <https://github.com/joonho112/bhfvar/issues>

## Author

JoonHo Lee <jlee296@ua.edu>
