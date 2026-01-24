# Prepare Data for BHF Model

Transforms survey data into the format required by the BHF Stan model.
This function handles index recoding, weight scaling, and computation of
design effect estimates needed for de-attenuation.

## Usage

``` r
prepare_bhf_data(
  data,
  outcome,
  domain,
  strata,
  psu,
  weights,
  population_shares = NULL,
  use_deattenuation = TRUE,
  prior_alpha_mean = NULL,
  prior_alpha_sd = 1.5
)
```

## Arguments

- data:

  A data frame containing the survey data.

- outcome:

  Character string. Name of the binary outcome variable (0/1).

- domain:

  Character string. Name of the domain/state variable.

- strata:

  Character string. Name of the stratification variable.

- psu:

  Character string. Name of the PSU (primary sampling unit) variable.

- weights:

  Character string. Name of the survey weight variable.

- population_shares:

  Optional numeric vector of length S (number of domains) containing
  population shares for each domain. Must sum to 1. If NULL, shares are
  estimated from the weighted data.

- use_deattenuation:

  Logical. If TRUE (default), computes and applies de-attenuation
  adjustment for finite-sample variance inflation.

- prior_alpha_mean:

  Numeric. Prior mean for the intercept on logit scale. Default is NULL,
  which estimates from data.

- prior_alpha_sd:

  Numeric. Prior SD for the intercept. Default is 1.5.

## Value

An object of class `bhf_data` containing:

- stan_data:

  List of data formatted for Stan

- mapping:

  List containing domain/strata/PSU label mappings

- domain_summary:

  Data frame with domain-level summary statistics

- input_info:

  List recording input column names and settings

## Details

This function performs several critical transformations:

- Index Recoding:

  All grouping variables are recoded to consecutive integers starting
  from 1 (required by Stan).

- Weight Scaling:

  Weights are scaled using Method D2 from Pfeffermann et al. (1998) to
  have effective sample size within each domain. This is critical for
  proper pseudo-likelihood estimation.

- Sampling Variance Estimation:

  For each domain, estimates the sampling variance of the proportion
  using the design effect.

- PSU Structure:

  Creates the nested PSU-within-stratum structure required by the Stan
  model.

## Weight Scaling (Method D2)

Weights are scaled so that for each domain s: \$\$w^\*\_i = w_i \times
\frac{n^{eff}\_s}{\sum\_{i \in s} w_i}\$\$ where \\n^{eff}\_s = (\sum
w_i)^2 / \sum w_i^2\\ is the effective sample size. This ensures the
pseudo-likelihood contributes appropriate information.

## Sampling Variance Estimation

The estimated sampling variance for domain s is: \$\$\hat{V}\_s =
\frac{deff_s \times \hat{p}\_s(1-\hat{p}\_s)}{n_s}\$\$ where \\deff_s\\
is the design effect. A default value of 1.5 is used when the design
effect cannot be reliably estimated.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example data
data(bhf_synthetic_data)

# Prepare data for Stan
prepared <- prepare_bhf_data(
  data = bhf_synthetic_data,
  outcome = "has_subsidy",
  domain = "state",
  strata = "stratum",
  psu = "psu",
  weights = "weight"
)

# Inspect the result
print(prepared)
summary(prepared)
} # }
```
