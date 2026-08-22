# Plot Domain Estimates

A caterpillar plot of the domain-level probabilities with their
posterior intervals, ordered by the posterior mean.

## Usage

``` r
bhf_plot_domains(fit, estimand = c("A", "B"), prob = 0.95, n_domains = NULL)
```

## Arguments

- fit:

  A `bhf_fit`.

- estimand:

  `"A"` or `"B"`.

- prob:

  Interval probability.

- n_domains:

  Optional cap on how many domains to draw, taking the extremes at each
  end. `NULL` draws all of them. For an odd cap, the extra domain is
  taken from the upper end of the ordered estimates.

## Value

A `ggplot` object.

## See also

[`domain_estimates()`](https://joonho112.github.io/bhfvar/reference/domain_estimates.md),
[`bhf_plot_shrinkage()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_shrinkage.md)

## Examples

``` r
if (FALSE) { # \dontrun{
bhf_plot_domains(fit, estimand = "B")
} # }
```
