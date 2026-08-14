# Resolve design-based sampling variances

Validates and aligns fixed, externally calculated domain sampling
variances, or records that de-attenuation is unavailable. Taylor
calculation is delegated to
[`estimate_taylor_vhat()`](https://joonho112.github.io/bhfvar/reference/estimate_taylor_vhat.md).

## Usage

``` r
resolve_sampling_variances(
  deattenuation = NULL,
  sampling_variances = NULL,
  sampling_variance_method = NULL,
  domain_labels,
  analysis_data = NULL,
  use_deattenuation = NULL
)
```
