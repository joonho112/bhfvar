# Plot Shrinkage of Domain Estimates

Compares the weighted raw proportion in each domain against the
model-based estimate and its posterior interval. Differences may reflect
partial pooling as well as standardization over the fitted survey-design
effects.

## Usage

``` r
bhf_plot_shrinkage(fit, estimand = c("A", "B"), prob = 0.95)
```

## Arguments

- fit:

  A `bhf_fit`.

- estimand:

  `"A"` or `"B"`.

- prob:

  Interval probability.

## Value

A `ggplot` object.

## See also

[`domain_estimates()`](https://joonho112.github.io/bhfvar/reference/domain_estimates.md),
[`bhf_plot_domains()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_domains.md)

## Examples

``` r
if (FALSE) { # \dontrun{
bhf_plot_shrinkage(fit)
} # }
```
