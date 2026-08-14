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
  weight_scaling = c("mean_one", "legacy_d2"),
  deattenuation = NULL,
  sampling_variances = NULL,
  sampling_variance_method = NULL,
  use_deattenuation = NULL,
  prior_alpha_mean = NULL,
  prior_alpha_sd = 0.5,
  sigma_state_prior = names(sigma_state_prior_catalog())
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

  Optional named positive numeric vector whose names exactly match the
  observed complete-case domains. Values are normalized to sum to one in
  canonical domain order. If NULL, shares are estimated from retained
  raw survey weights.

- weight_scaling:

  Likelihood-weight scaling method. The default `"mean_one"` scales raw
  weights globally to sum to the retained sample size. `"legacy_d2"` is
  a warned sensitivity/reproduction path.

- deattenuation:

  Character mode: `"taylor"` (default), `"supplied"`, or `"none"`.

- sampling_variances:

  Named nonnegative domain sampling variances for
  `deattenuation = "supplied"`.

- sampling_variance_method:

  Provenance label for supplied variances: `"external_taylor"`,
  `"external_replicate"`, or `"external_other"`.

- use_deattenuation:

  Deprecated logical compatibility argument. Use `deattenuation`
  instead.

- prior_alpha_mean:

  Numeric. Prior mean for the intercept on logit scale. Default is NULL,
  which estimates from data.

- prior_alpha_sd:

  Numeric. Prior SD for the intercept. Default is 0.5.

- sigma_state_prior:

  Article prior-sensitivity selector for the state random-effect SD.
  Exactly one of `"half_t3_2.5"` (baseline), `"half_normal_1"`,
  `"half_cauchy_2.5"`, or `"half_t3_5"`. Stratum and PSU SD priors
  remain half-t(3, 0, 2.5).

## Value

An object of class `bhf_data` containing:

- schema_version:

  Prepared-data schema version.

- contract_id:

  Stable prepared-data contract identifier.

- stan_data:

  List of data formatted for Stan

- mapping:

  List containing domain/strata/PSU label mappings

- row_provenance:

  Original/retained/dropped row mapping and ledger.

- provenance:

  Consolidated row, weight, share, vhat, and prior metadata.

- analysis_data:

  Compact retained analysis frame with raw and likelihood weights.

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

  Positive raw weights are retained and globally normalized once in R to
  have mean one.

- Sampling Variance Estimation:

  For each domain, estimates the sampling variance of the proportion
  using the design effect.

- PSU Structure:

  Creates the nested PSU-within-stratum structure required by the Stan
  model.

## Weight Scaling

The default likelihood weight is \$\$w^\*\_i = w_i \times
\frac{N}{\sum_i w_i},\$\$ so the scaled weights sum to the retained
sample size N while preserving all raw-weight ratios. Legacy D2 scaling
is available only through the explicit warned
`weight_scaling = "legacy_d2"` sensitivity path; it receives a final
global sum-N normalization in R.

## Sampling-variance contract

The default `deattenuation = "taylor"` uses the retained raw survey
weights and the one-stage stratified PSU design. Singleton strata fail
explicitly; no cap, floor, or default design effect is substituted.
`deattenuation = "supplied"` treats named domain variances as fixed
external inputs. Their estimation uncertainty is not propagated.

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
