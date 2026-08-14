# Extract Raw or Pseudo Pointwise Log-Likelihood Contributions

Aggregation only sums stored observation contributions. It does not
refit the model and does not establish ordinary observation- or
cluster-level LOO.

## Usage

``` r
log_lik(
  fit,
  kind = c("pseudo", "raw"),
  aggregate = c("observation", "psu", "stratum")
)
```

## Arguments

- fit:

  A versioned `bhf_fit` object.

- kind:

  `"pseudo"` (default) or `"raw"`.

- aggregate:

  `"observation"`, `"psu"`, or `"stratum"`.

## Value

A draw-by-unit `bhf_log_lik` matrix with explicit `kind`, `aggregate`,
and interpretation metadata.

## Details

The default pseudo contributions target the fitted pseudo-likelihood.
Aggregation is arithmetic only; this function does not establish
ordinary observation- or cluster-level LOO validity.
