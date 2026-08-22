# Sensitivity of Estimand A\* to the Sampling Variances

Recomputes Estimand A\* with the supplied sampling variances multiplied
by each of several factors. Because A\* is a deterministic function of
Estimand A's posterior draws and the fixed \\\hat V_s\\, this needs no
new sampling.

## Usage

``` r
bhf_astar_sensitivity(fit, scale = c(0.5, 0.75, 1, 1.25, 1.5), prob = 0.95)
```

## Arguments

- fit:

  A `bhf_fit` fitted with de-attenuation enabled.

- scale:

  Numeric multipliers applied to the supplied \\\hat V_s\\.

- prob:

  Interval probability.

## Value

A data frame of class `bhf_astar_sensitivity`, one row per scale, with
the posterior mean, interval, and boundary diagnostics for A\*'s
between-domain variance and proportion.

## Details

The article treats the \\\hat V_s\\ as known and does not propagate
their estimation uncertainty (article Section 2.4). This function does
not change that assumption. It reports how sensitive the conclusion is
to it. A large posterior share at the zero boundary indicates a
boundary-dominated, nonregular A\* summary; it does not make A\*
mathematically undefined.

## See also

[`variance_decomposition()`](https://joonho112.github.io/bhfvar/reference/variance_decomposition.md),
[`bhf_plot_astar_sensitivity()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_astar_sensitivity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
bhf_astar_sensitivity(fit)
bhf_astar_sensitivity(fit, scale = c(0.5, 1, 2))
} # }
```
