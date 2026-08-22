# Plot the Probability-Scale Variance Decomposition

Shows the between-domain, within-domain, and total variance for
Estimands A, A\*, and B side by side, so the effect of the design (A
versus B) and of the de-attenuation correction (A versus A\*) are
visible at once.

## Usage

``` r
bhf_plot_variance(x, components = c("between", "within", "total"), prob = 0.95)
```

## Arguments

- x:

  A `bhf_fit`, or the result of
  [`variance_decomposition()`](https://joonho112.github.io/bhfvar/reference/variance_decomposition.md).

- components:

  Which variance components to show.

- prob:

  Interval probability, used when `x` is a `bhf_fit`.

## Value

A `ggplot` object.

## See also

[`variance_decomposition()`](https://joonho112.github.io/bhfvar/reference/variance_decomposition.md),
[`bhf_plot_icc()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_icc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- bhf_fit(prepared, model = model)
bhf_plot_variance(fit)
} # }
```
