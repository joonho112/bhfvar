# Generate a Small Article-Aligned Synthetic Survey

Internal deterministic generator for scientific validation. It
implements the published low/high variance profiles and
informative-weight tilt without running the article's large Monte Carlo
study.

## Usage

``` r
generate_article_synthetic(
  profile = c("low", "high"),
  rho = c(0, 0.5),
  n_states = 4L,
  n_strata = 3L,
  psus_per_stratum = 2L,
  observations_per_cell = 2L,
  alpha = -1.5,
  population_shares = NULL,
  vhat_state = NULL,
  effect_seed = 7106L,
  outcome_seed = 7107L
)
```
