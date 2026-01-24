# Scale Weights Using Method D2

Internal function to scale survey weights using the D2 method
(Pfeffermann et al., 1998).

## Usage

``` r
scale_weights_d2(weight, state_id, domain_summary)
```

## Arguments

- weight:

  Original survey weights.

- state_id:

  State/domain IDs.

- domain_summary:

  Domain summary from compute_domain_summary().

## Value

Scaled weights.
