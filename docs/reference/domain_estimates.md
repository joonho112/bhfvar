# Extract Domain-Level Estimand A or B Probabilities

Extract Domain-Level Estimand A or B Probabilities

## Usage

``` r
domain_estimates(fit, estimand = c("A", "B"), prob = 0.95, type = NULL)
```

## Arguments

- fit:

  A versioned `bhf_fit` object.

- estimand:

  Either `"A"` or `"B"`.

- prob:

  Credible-interval probability.

- type:

  Deprecated transition argument. `"conditional"` maps to A with a
  classed warning; `"marginal"` fails because it is not Estimand B.

## Value

A `bhf_domain_estimates` data frame in canonical domain-ID order, with
label/ID, population share, sample size, posterior summaries, dynamic
interval bounds, and estimand metadata.
