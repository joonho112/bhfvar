# Compute MCMC Diagnostics

Internal function to compute diagnostic summaries from a stanfit object.

## Usage

``` r
compute_diagnostics(
  stanfit,
  max_treedepth = 12L,
  pars = c("alpha", "sigma_state", "sigma_stratum", "sigma_psu")
)
```

## Arguments

- stanfit:

  A stanfit object.

- max_treedepth:

  Configured maximum tree depth used for this fit.

- pars:

  Structural parameter names to diagnose. Generated quantities are
  intentionally excluded.

## Value

List of diagnostic summaries.
