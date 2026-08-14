# Introduction to bhfvar

## Purpose

`bhfvar` is an experimental implementation of a Bayesian
pseudo-likelihood framework for decomposing binary-outcome variation
across post-hoc domains in a complex survey. It separates a substantive
domain effect from crossed stratum and PSU-within-stratum design
effects.

The package addresses two distinct questions:

1.  What probability-scale variation is induced by domain effects alone
    (A)?
2.  What variation is present after retaining the observed mixture of
    design effects within each domain (B)?

A\* subtracts a fixed, population-share-weighted sampling-variance
correction from A’s between-domain component, truncated at zero.

## Model structure

For observation $`i`$,

``` math
\operatorname{logit}(p_i)=\alpha+u_{s[i]}+v_{h[i]}+b_{j[i]},
```

where domain and stratum effects are crossed and each flat PSU
identifier is nested in exactly one stratum. State effects are centered
using population shares, stratum effects by their simple mean, and PSU
effects within stratum.

``` text
observations
 ├─ domain/state (crossed)
 └─ stratum
     └─ PSU (nested)
```

## Supported workflow

``` r

library(bhfvar)
data(bhf_synthetic_data)

prepared <- prepare_bhf_data(
  bhf_synthetic_data,
  outcome = "has_subsidy",
  domain = "state",
  strata = "stratum",
  psu = "psu",
  weights = "weight",
  deattenuation = "taylor"
)
model <- compile_bhf_model()
fit <- bhf_fit(prepared, model = model, chains = 4, iter = 2000,
               seed = 1234)

variance_decomposition(fit)
domain_estimates(fit, estimand = "A")
domain_estimates(fit, estimand = "B")
```

The example is intentionally not evaluated while building the vignette
because Stan compilation and MCMC require a configured toolchain. The
data-only path is evaluated in the quick-start vignette and dedicated
tests cover public result contracts.

## Scope

The implementation is verified against the model specified in Appendix A
of the accompanying article: a frozen reference oracle, algebraic
property tests, and centering-invariant tests all pass.

Interval calibration is a separate question and has not been
established. Posterior intervals here are pseudo-posterior credible
intervals; their frequentist coverage is not guaranteed, and a variance
adjustment is not implemented in this version. The article’s
restricted-data application is not reproduced by this package.

A successful compile, fit, or extraction shows that the workflow ran. It
says nothing about interval calibration. See
[`vignette("scientific-limitations")`](https://joonho112.github.io/bhfvar/articles/scientific-limitations.md).

## Scope boundaries

- A\* conditions on supplied or Taylor-estimated sampling variances;
  uncertainty in those estimates is not propagated.
- Pseudo-posterior intervals require design-aware interpretation.
- [`log_lik()`](https://joonho112.github.io/bhfvar/reference/log_lik.md)
  does not establish ordinary observation- or cluster-level LOO.
- No plotting API or posterior-predictive helper is supported.
- The bundled legacy synthetic data are illustrative, not restricted
  NSECE data and not the confirmatory recovery fixture.

See
[`vignette("methodology", package = "bhfvar")`](https://joonho112.github.io/bhfvar/articles/methodology.md),
[`vignette("dual-estimands", package = "bhfvar")`](https://joonho112.github.io/bhfvar/articles/dual-estimands.md),
and
[`vignette("diagnostics", package = "bhfvar")`](https://joonho112.github.io/bhfvar/articles/diagnostics.md)
before substantive use.
