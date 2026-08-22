# Plot Latent-Scale Intraclass Correlations

A latent-scale diagnostic formed from the model's pre-centering
random-effect scale parameters and the logistic residual variance. These
shares are not empirical variances of the centered realized effects and
do not by themselves determine whether a random-effect term is
warranted.

## Usage

``` r
bhf_plot_icc(x, prob = 0.95)
```

## Arguments

- x:

  A `bhf_fit`, or the result of
  [`variance_decomposition()`](https://joonho112.github.io/bhfvar/reference/variance_decomposition.md).

- prob:

  Interval probability, used when `x` is a `bhf_fit`.

## Value

A `ggplot` object.

## See also

[`variance_decomposition()`](https://joonho112.github.io/bhfvar/reference/variance_decomposition.md),
[`bhf_plot_variance()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_variance.md)

## Examples

``` r
if (FALSE) { # \dontrun{
bhf_plot_icc(fit)
} # }
```
