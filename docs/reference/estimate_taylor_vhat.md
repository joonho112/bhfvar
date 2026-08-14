# Estimate Domain Sampling Variances by Taylor Linearization

Uses a one-stage, stratified PSU survey design and raw survey weights to
estimate the design-based variance of each domain Hajek proportion. The
function deliberately fails for singleton strata and does not apply
caps, floors, or heuristic fallbacks.

## Usage

``` r
estimate_taylor_vhat(analysis_data, domain_labels)
```
