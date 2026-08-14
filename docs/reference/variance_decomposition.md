# Extract Article-Aligned Variance Decomposition Results

Summarizes separate latent diagnostics and probability-scale Estimands
A, A\*, and B. No legacy 0.3.0 quantity is silently reinterpreted.

## Usage

``` r
variance_decomposition(fit, prob = 0.95, print = FALSE)
```

## Arguments

- fit:

  A versioned `bhf_fit` object with schema `0.5.0`.

- prob:

  Credible-interval probability.

- print:

  Retained transition argument. If `TRUE`, prints the tidy table.

## Value

A `bhf_variance_decomposition` list with `latent`, `A`, `A_star`, `B`,
`gaps`, and `summary_table` components. `A_star$available` is `FALSE`
when de-attenuation was disabled.

## Details

Estimands A, A\*, and B are probability-scale decompositions;
latent-scale SDs, variances, and ICCs are reported separately. A\*
treats sampling variances as fixed inputs and does not propagate their
estimation uncertainty. Intervals are pseudo-posterior credible
intervals; their frequentist coverage has not been established.
