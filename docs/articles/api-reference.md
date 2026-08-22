# Public API and Object Contracts

## Exported functions

| Function | Current signature |
|----|----|
| `compile_bhf_model` | `(verbose = TRUE, auto_write = TRUE)` |
| `get_stan_file_path` | `()` |
| `prepare_bhf_data` | `(data, outcome, domain, strata, psu, weights, population_shares = NULL, weight_scaling = c("mean_one", "legacy_d2"), deattenuation = NULL, sampling_variances = NULL, sampling_variance_method = NULL, use_deattenuation = NULL, prior_alpha_mean = NULL, prior_alpha_sd = 0.5, sigma_state_prior = names(sigma_state_prior_catalog()))` |
| `bhf_fit` | `(data, model = NULL, chains = 4, iter = 2000, warmup = floor(iter/2), seed = 1234, cores = NULL, adapt_delta = 0.95, max_treedepth = 12, refresh = 200, ...)` |
| `variance_decomposition` | `(fit, prob = 0.95, print = FALSE)` |
| `domain_estimates` | `(fit, estimand = c("A", "B"), prob = 0.95, type = NULL)` |
| `overall_estimate` | `(fit, estimand = c("A", "B"), prob = 0.95)` |
| `log_lik` | `(fit, kind = c("pseudo", "raw"), aggregate = c("observation", "psu", "stratum"))` |
| `bhf_plot_variance` | `(x, components = c("between", "within", "total"), prob = 0.95)` |
| `bhf_plot_icc` | `(x, prob = 0.95)` |
| `bhf_plot_domains` | `(fit, estimand = c("A", "B"), prob = 0.95, n_domains = NULL)` |
| `bhf_plot_shrinkage` | `(fit, estimand = c("A", "B"), prob = 0.95)` |
| `bhf_astar_sensitivity` | `(fit, scale = c(0.5, 0.75, 1, 1.25, 1.5), prob = 0.95)` |
| `bhf_plot_astar_sensitivity` | `(x, what = c("proportion", "between"))` |
| `calc_eff_n` | `(weights)` |

Function help pages are authoritative for argument-level details.

## Versioned objects

[`prepare_bhf_data()`](https://joonho112.github.io/bhfvar/reference/prepare_bhf_data.md)
returns a `bhf_data` object with schema, contract ID, Stan data,
canonical mappings, row provenance, analysis data, domain summary, and
weight/population-share/sampling-variance/prior provenance.

[`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md)
returns a `bhf_fit` with schema `0.5.0`, model hash, `stanfit`, the
prepared data, diagnostics, sampling call, and consolidated provenance.
Legacy objects fail closed or are identified by internal migration
checks; they are not silently relabeled as 0.5.0.

## Extractor returns

- [`variance_decomposition()`](https://joonho112.github.io/bhfvar/reference/variance_decomposition.md)
  returns `latent`, `A`, `A_star`, `B`, `gaps`, and a tidy
  `summary_table`.
- [`domain_estimates()`](https://joonho112.github.io/bhfvar/reference/domain_estimates.md)
  returns one row per canonical domain with dynamic interval bounds and
  population shares.
- [`overall_estimate()`](https://joonho112.github.io/bhfvar/reference/overall_estimate.md)
  returns a population-share-weighted probability summary.
- [`log_lik()`](https://joonho112.github.io/bhfvar/reference/log_lik.md)
  returns a draw-by-unit matrix with explicit scope metadata.

All interval labels are derived from `prob`; code must not assume fixed
column names such as `q2.5` or `q97.5`.

## Deprecations

- `weight_scaling = "legacy_d2"` warns and is retained only for
  transition and sensitivity comparison.
- `use_deattenuation` is replaced by `deattenuation`.
- `domain_estimates(type = "conditional")` warns and maps to A.
- `domain_estimates(type = "marginal")` errors; it is not B.
- 0.3.0 result fields are not silently redefined.

## Remaining unsupported surfaces

There is no reliability statistic, sandwich adjustment, ordinary LOO
helper, or posterior-predictive-check helper. The exported plot helpers
below return `ggplot` objects when the suggested `ggplot2` package is
installed.

## Plots and diagnostics

``` r

# Probability-scale decomposition across A, A*, and B
bhf_plot_variance(fit)

# Latent shares formed from the pre-centering random-effect scale parameters
bhf_plot_icc(fit)

# Domain estimates, ordered, with intervals; n_domains keeps the extremes
bhf_plot_domains(fit, estimand = "B", n_domains = 20)

# Model-based estimate against the weighted raw proportion
bhf_plot_shrinkage(fit)
```

All four need `ggplot2`, a suggested dependency. Each returns a `ggplot`
object, so you can add layers or themes as usual.

Estimand A\* subtracts a fixed sampling-variance correction and does not
propagate the uncertainty in those variances. To see what that
assumption is worth:

``` r

s <- bhf_astar_sensitivity(fit, scale = c(0.5, 0.75, 1, 1.25, 1.5))
s
bhf_plot_astar_sensitivity(s)
```

The `at zero` column is the posterior share of draws truncated at the
zero boundary. A large share indicates a boundary-dominated, nonregular
A\* summary; it does not make A\* mathematically undefined.

## Scope

The API contracts and frozen-oracle tests pass within their scope: the
functions return what they document, and the computed quantities match
the article’s definitions. Repeated-sample coverage has been measured
only in a reduced balanced synthetic design and is not guaranteed more
generally — see
[`vignette("scientific-limitations")`](https://joonho112.github.io/bhfvar/articles/scientific-limitations.md).
