# Resolve Population Shares for the Observed Domain Universe

Internal helper that validates externally supplied population shares or
estimates them from retained raw survey weights. Returned shares are
named, aligned to `domain_labels`, strictly positive, and normalized to
sum to one.

## Usage

``` r
resolve_population_shares(
  population_shares,
  domain_labels,
  state_id,
  raw_weights
)
```

## Arguments

- population_shares:

  NULL or a named numeric vector.

- domain_labels:

  Observed domain labels in canonical state-ID order.

- state_id:

  Integer state IDs aligned with `raw_weights`.

- raw_weights:

  Positive finite raw survey weights.

## Value

An internal `bhf_population_shares` object.
