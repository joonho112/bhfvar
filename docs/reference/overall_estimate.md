# Extract Overall Population Estimate

Extracts the posterior distribution of the overall (population-weighted)
probability from a fitted BHF model.

## Usage

``` r
overall_estimate(fit, prob = 0.95)
```

## Arguments

- fit:

  An object of class `bhf_fit` from
  [`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md).

- prob:

  Numeric. Probability for credible interval. Default is 0.95.

## Value

A list with components:

- mean:

  Posterior mean

- sd:

  Posterior standard deviation

- q025:

  Lower credible interval bound

- q500:

  Posterior median

- q975:

  Upper credible interval bound
