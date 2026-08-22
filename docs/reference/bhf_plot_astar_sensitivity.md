# Plot the A\* Sampling-Variance Sensitivity

Plot the A\* Sampling-Variance Sensitivity

## Usage

``` r
bhf_plot_astar_sensitivity(x, what = c("proportion", "between"))
```

## Arguments

- x:

  The result of
  [`bhf_astar_sensitivity()`](https://joonho112.github.io/bhfvar/reference/bhf_astar_sensitivity.md).

- what:

  Which quantity to plot.

## Value

A `ggplot` object.

## See also

[`bhf_astar_sensitivity()`](https://joonho112.github.io/bhfvar/reference/bhf_astar_sensitivity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
bhf_plot_astar_sensitivity(bhf_astar_sensitivity(fit))
} # }
```
