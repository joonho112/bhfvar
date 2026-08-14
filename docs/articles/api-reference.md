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

## Unsupported surfaces

There is no exported plotting method, reliability statistic, sandwich
adjustment, ordinary LOO helper, or posterior-predictive-check helper.
Use tidy outputs for local visualization without attributing a
package-level plot API.

## Scope

The API contracts and frozen-oracle tests pass within their scope: the
functions return what they document, and the computed quantities match
the article’s definitions. Interval calibration is a separate matter and
has not been established — see
[`vignette("scientific-limitations")`](https://joonho112.github.io/bhfvar/articles/scientific-limitations.md).
