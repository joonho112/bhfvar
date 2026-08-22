# Quick Start

## A checked data-only path

This chunk is evaluated during vignette build. De-attenuation is
disabled here so the example remains lightweight and does not turn a
legacy illustrative data object into evidence about sampling-variance
recovery.

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
  weight_scaling = "mean_one",
  deattenuation = "none"
)

list(
  schema = prepared$schema_version,
  observations = prepared$stan_data$N,
  domains = prepared$stan_data$S,
  strata = prepared$stan_data$H,
  nested_psus = prepared$stan_data$J
)
#> $schema
#> [1] "0.5.0"
#> 
#> $observations
#> [1] 1598
#> 
#> $domains
#> [1] 50
#> 
#> $strata
#> [1] 27
#> 
#> $nested_psus
#> [1] 356
sum(prepared$stan_data$w_lik)
#> [1] 1598
```

The final line equals the retained sample size: default likelihood
weights are globally normalized to mean one.

## Compile and fit locally

The following is not evaluated during documentation builds.

``` r

model <- compile_bhf_model()
fit <- bhf_fit(
  prepared,
  model = model,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 1234,
  adapt_delta = 0.95,
  max_treedepth = 12
)
print(fit)
fit$diagnostics
```

`fit$diagnostics` reports divergences, configured/observed tree depth,
maximum R-hat, and legacy rstan `n_eff`. It does not label `n_eff` as
bulk or tail ESS.

## Extract explicit quantities

``` r

vd <- variance_decomposition(fit, prob = 0.95)
vd$A$summary
vd$A_star$available
vd$A_star$provenance
vd$B$summary
vd$gaps$table

domain_estimates(fit, estimand = "A", prob = 0.95)
domain_estimates(fit, estimand = "B", prob = 0.95)
overall_estimate(fit, estimand = "B", prob = 0.95)
```

If `deattenuation = "none"`, A\* remains in the result schema but
`vd$A_star$available` is `FALSE`. Re-prepare with
`deattenuation = "taylor"` or `"supplied"` to compute A\* under the
documented fixed-input limitation.

## Interpretation checkpoint

A workflow that runs cleanly tells you the pipeline worked; it tells you
nothing by itself about interval calibration. The intervals returned
here are pseudo-posterior credible intervals. Repeated-sample coverage
has been measured only in a reduced balanced synthetic design and is not
guaranteed more generally. See
[`vignette("scientific-limitations")`](https://joonho112.github.io/bhfvar/articles/scientific-limitations.md)
before reporting them.
