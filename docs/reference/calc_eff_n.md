# Calculate Effective Sample Size

Calculates the effective sample size given survey weights.

## Usage

``` r
calc_eff_n(weights)
```

## Arguments

- weights:

  Numeric vector of survey weights.

## Value

Effective sample size (scalar).

## Examples

``` r
w <- c(1, 2, 1.5, 3, 2)
calc_eff_n(w)  # Should be less than length(w)
#> [1] 4.45679
```
