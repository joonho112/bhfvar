# Extract an Overall Estimand A or B Population Probability

Extract an Overall Estimand A or B Population Probability

## Usage

``` r
overall_estimate(fit, estimand = c("A", "B"), prob = 0.95)
```

## Arguments

- fit:

  A versioned `bhf_fit` object.

- estimand:

  Either `"A"` or `"B"`.

- prob:

  Credible-interval probability.

## Value

A `bhf_overall_estimate` list containing the estimand, interval
contract, posterior summary, and population-share provenance.
