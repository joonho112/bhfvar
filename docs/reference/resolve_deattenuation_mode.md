# Resolve the de-attenuation mode

Internal compatibility helper. A `NULL` new-style argument means that
the caller did not explicitly select a mode, so the article-aligned
Taylor mode is the default. The deprecated logical argument is retained
for the 0.4.x transition only.

## Usage

``` r
resolve_deattenuation_mode(
  deattenuation = NULL,
  sampling_variances = NULL,
  use_deattenuation = NULL
)
```
