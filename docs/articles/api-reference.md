# API Reference: Complete Function Documentation

## Overview

This vignette provides comprehensive documentation of all functions
exported by the bhfvar package. Functions are organized by workflow
stage:

1.  **Model Compilation**
2.  **Data Preparation**
3.  **Model Fitting**
4.  **Result Extraction**
5.  **Utility Functions**

## 1. Model Compilation

### compile_bhf_model

Compiles the Stan model for the Bayesian Hybrid Framework.

**Usage**

``` r
compile_bhf_model(verbose = TRUE, auto_write = TRUE)
```

**Arguments**

| Argument     | Type    | Default | Description                  |
|--------------|---------|---------|------------------------------|
| `verbose`    | logical | TRUE    | Print compilation progress   |
| `auto_write` | logical | TRUE    | Cache compiled model to disk |

**Value**

An object of class `stanmodel` that can be passed to
[`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md).

**Details**

This function should be called once at the beginning of each R session.
The first compilation takes 30–60 seconds; subsequent calls use the
cached version (if `auto_write = TRUE`).

**Examples**

``` r
library(bhfvar)

# Compile with default settings
model <- compile_bhf_model()

# Silent compilation
model <- compile_bhf_model(verbose = FALSE)
```

### get_stan_file_path

Returns the file path to the Stan model file.

**Usage**

``` r
get_stan_file_path()
```

**Value**

A character string with the full path to `bhf_hybrid.stan`.

**Examples**

``` r
# Get the path
stan_path <- get_stan_file_path()
cat(stan_path)

# Read the Stan code
stan_code <- readLines(stan_path)
head(stan_code, 20)
```

## 2. Data Preparation

### prepare_bhf_data

Transforms survey data into the format required by the Stan model.

**Usage**

``` r
prepare_bhf_data(
  data,
  outcome,
  domain,
  strata,
  psu,
  weights = NULL,
  use_deattenuation = TRUE
)
```

**Arguments**

| Argument | Type | Default | Description |
|----|----|----|----|
| `data` | data.frame | — | Input data |
| `outcome` | character | — | Name of binary outcome column |
| `domain` | character | — | Name of domain indicator column |
| `strata` | character | — | Name of stratum indicator column |
| `psu` | character | — | Name of PSU indicator column |
| `weights` | character | NULL | Name of weight column (NULL = equal weights) |
| `use_deattenuation` | logical | TRUE | Compute sampling variances for de-attenuation |

**Value**

A list of class `bhf_data` containing:

- `stan_data`: List formatted for Stan
- `domain_map`: Mapping from domain names to indices
- `stratum_map`: Mapping from stratum names to indices
- `psu_map`: Mapping from PSU names to indices
- `sampling_variances`: Vector of $`\hat{V}_s`$ (if
  `use_deattenuation = TRUE`)
- `original_data`: Copy of input data

**Details**

The function:

1.  Validates input data (checks for missing values, proper types)
2.  Creates integer indices for domains, strata, and PSUs
3.  Computes population shares for each domain
4.  Scales weights using Method 2 (sum to sample size)
5.  Estimates design-based sampling variances via Taylor linearization
6.  Formats everything for Stan

**Examples**

``` r
library(bhfvar)
data(bhf_synthetic_data)

# Basic preparation
prepared <- prepare_bhf_data(
  data = bhf_synthetic_data,
  outcome = "has_subsidy",
  domain = "state",
  strata = "stratum",
  psu = "psu",
  weights = "weight"
)

print(prepared)

# Check dimensions
cat("N:", prepared$stan_data$N, "\n")
cat("S:", prepared$stan_data$S, "\n")
cat("H:", prepared$stan_data$H, "\n")
cat("J:", prepared$stan_data$J, "\n")

# Disable de-attenuation
prepared_no_deatten <- prepare_bhf_data(
  bhf_synthetic_data, "has_subsidy", "state", "stratum", "psu", "weight",
  use_deattenuation = FALSE
)
```

## 3. Model Fitting

### bhf_fit

Fits the Bayesian Hybrid Framework model via Stan.

**Usage**

``` r
bhf_fit(
  data,
  model,
  chains = 4,
  iter = 2000,
  warmup = floor(iter/2),
  thin = 1,
  seed = NULL,
  control = list(),
  ...
)
```

**Arguments**

| Argument | Type | Default | Description |
|----|----|----|----|
| `data` | bhf_data | — | Output from [`prepare_bhf_data()`](https://joonho112.github.io/bhfvar/reference/prepare_bhf_data.md) |
| `model` | stanmodel | — | Output from [`compile_bhf_model()`](https://joonho112.github.io/bhfvar/reference/compile_bhf_model.md) |
| `chains` | integer | 4 | Number of MCMC chains |
| `iter` | integer | 2000 | Total iterations per chain |
| `warmup` | integer | iter/2 | Warmup iterations (discarded) |
| `thin` | integer | 1 | Thinning interval |
| `seed` | integer | NULL | Random seed for reproducibility |
| `control` | list | list() | Control parameters (e.g., `adapt_delta`) |
| `...` | — | — | Additional arguments passed to [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) |

**Value**

A list of class `bhf_fit` containing:

- `stanfit`: The `stanfit` object from rstan
- `data`: The prepared data used for fitting
- `diagnostics`: Convergence diagnostics (divergences, Rhat, ESS)
- `settings`: Sampling settings used

**Details**

The function runs Stan’s NUTS sampler in parallel (using `mc.cores`
available cores). Convergence diagnostics are automatically computed and
warnings issued if problems are detected.

**Examples**

``` r
library(bhfvar)
library(rstan)

# Setup
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Compile and prepare
model <- compile_bhf_model()
data(bhf_synthetic_data)
prepared <- prepare_bhf_data(
  bhf_synthetic_data, "has_subsidy", "state", "stratum", "psu", "weight"
)

# Basic fit
fit <- bhf_fit(prepared, model = model)

# With more iterations
fit <- bhf_fit(prepared, model = model, chains = 4, iter = 4000)

# With control parameters for divergences
fit <- bhf_fit(prepared, model = model,
               control = list(adapt_delta = 0.95, max_treedepth = 12))

# Reproducible
fit <- bhf_fit(prepared, model = model, seed = 42)
```

## 4. Result Extraction

### variance_decomposition

Extracts variance decomposition across all three estimands.

**Usage**

``` r
variance_decomposition(fit, prob = 0.95, print = TRUE)
```

**Arguments**

| Argument | Type | Default | Description |
|----|----|----|----|
| `fit` | bhf_fit | — | Output from [`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md) |
| `prob` | numeric | 0.95 | Probability mass for credible intervals |
| `print` | logical | TRUE | Print formatted summary |

**Value**

A list with components:

- `logit`: Estimand A results (logit scale)
  - `var_between_mean`, `var_between_q025`, `var_between_q975`
  - `var_within_mean`, `var_within_q025`, `var_within_q975`
  - `icc_mean`, `icc_q025`, `icc_q975`
- `prob`: Estimand B results (probability scale)
  - Same structure as `logit`
- `deatten`: Estimand A\* results (de-attenuated)
  - Same structure as `logit`
- `summary_table`: Formatted data frame of results

**Examples**

``` r
# Basic extraction
vd <- variance_decomposition(fit)

# Suppress printing
vd <- variance_decomposition(fit, print = FALSE)

# 90% credible intervals
vd <- variance_decomposition(fit, prob = 0.90)

# Access specific values
vd$prob$icc_mean      # Estimand B ICC
vd$deatten$icc_mean   # Estimand A* ICC
```

### domain_estimates

Extracts domain-specific probability estimates.

**Usage**

``` r
domain_estimates(fit, type = "marginal", prob = 0.95)
```

**Arguments**

| Argument | Type | Default | Description |
|----|----|----|----|
| `fit` | bhf_fit | — | Output from [`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md) |
| `type` | character | “marginal” | “marginal” or “conditional” |
| `prob` | numeric | 0.95 | Probability mass for credible intervals |

**Value**

A data frame with columns:

| Column        | Description                   |
|---------------|-------------------------------|
| `domain`      | Domain name                   |
| `domain_id`   | Domain index                  |
| `mean`        | Posterior mean probability    |
| `sd`          | Posterior standard deviation  |
| `q025`        | Lower credible interval bound |
| `q500`        | Median                        |
| `q975`        | Upper credible interval bound |
| `pop_share`   | Population share $`\pi_s`$    |
| `reliability` | Shrinkage factor $`R_s`$      |

**Details**

- **Marginal** estimates integrate out within-domain (PSU) variation
  using the Zeger adjustment. These are appropriate for population
  inference.

- **Conditional** estimates are $`\text{logit}^{-1}(\alpha + u_s)`$
  without integration. These represent the domain-specific propensity.

**Examples**

``` r
# Marginal estimates (default)
estimates <- domain_estimates(fit)
head(estimates)

# Conditional estimates
estimates_cond <- domain_estimates(fit, type = "conditional")

# Sort by estimate
estimates <- estimates[order(estimates$mean), ]

# Filter high-reliability domains
reliable <- estimates[estimates$reliability > 0.5, ]
```

### overall_estimate

Extracts the population-weighted overall probability.

**Usage**

``` r
overall_estimate(fit, prob = 0.95)
```

**Arguments**

| Argument | Type | Default | Description |
|----|----|----|----|
| `fit` | bhf_fit | — | Output from [`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md) |
| `prob` | numeric | 0.95 | Probability mass for credible intervals |

**Value**

A list with:

- `mean`: Posterior mean of overall probability
- `sd`: Posterior standard deviation
- `q025`, `q975`: Credible interval bounds
- `samples`: Full posterior samples (if needed for further analysis)

**Examples**

``` r
overall <- overall_estimate(fit)

cat("Overall:", round(overall$mean, 3),
    "(95% CI:", round(overall$q025, 3), "-", round(overall$q975, 3), ")\n")
```

### log_lik

Extracts the log-likelihood matrix for LOO-CV.

**Usage**

``` r
log_lik(fit)
```

**Arguments**

| Argument | Type | Default | Description |
|----|----|----|----|
| `fit` | bhf_fit | — | Output from [`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md) |

**Value**

A matrix of dimension (n_iterations × N) containing pointwise
log-likelihood values.

**Examples**

``` r
library(loo)

# Extract log-likelihood
ll <- log_lik(fit)

# Compute LOO
loo_result <- loo(ll)
print(loo_result)

# Compare models
loo_compare(loo_result, loo_result2)
```

## 5. Utility Functions

### validate_input_data

Validates input data before preparation.

**Usage**

``` r
validate_input_data(data, outcome, domain, strata, psu, weights = NULL)
```

**Arguments**

Same as
[`prepare_bhf_data()`](https://joonho112.github.io/bhfvar/reference/prepare_bhf_data.md).

**Value**

Returns `TRUE` invisibly if validation passes. Throws an error with
informative message if validation fails.

**Examples**

``` r
# Check data before preparation
validate_input_data(
  bhf_synthetic_data, "has_subsidy", "state", "stratum", "psu", "weight"
)
```

### validate_stan_data

Validates Stan data list before fitting.

**Usage**

``` r
validate_stan_data(stan_data)
```

**Arguments**

| Argument | Type | Default | Description |
|----|----|----|----|
| `stan_data` | list | — | Stan data list from [`prepare_bhf_data()`](https://joonho112.github.io/bhfvar/reference/prepare_bhf_data.md) |

**Value**

Returns `TRUE` invisibly if validation passes.

**Examples**

``` r
prepared <- prepare_bhf_data(
  bhf_synthetic_data, "has_subsidy", "state", "stratum", "psu", "weight"
)

validate_stan_data(prepared$stan_data)
```

## 6. Print and Summary Methods

### print.bhf_data

Print method for prepared data objects.

**Usage**

``` r
print(x, ...)
```

**Examples**

``` r
prepared <- prepare_bhf_data(
  bhf_synthetic_data, "has_subsidy", "state", "stratum", "psu", "weight"
)
print(prepared)
```

### print.bhf_fit

Print method for fitted model objects.

**Usage**

``` r
print(x, ...)
```

**Examples**

``` r
print(fit)
```

### summary.bhf_fit

Summary method providing variance decomposition and diagnostics.

**Usage**

``` r
summary(object, ...)
```

**Examples**

``` r
summary(fit)
```

## 7. Example Data

### bhf_synthetic_data

Synthetic dataset mimicking the NSECE structure.

**Description**

A data frame with 1,598 observations of 5 variables:

| Variable      | Type    | Description                    |
|---------------|---------|--------------------------------|
| `state`       | factor  | State identifier (50 states)   |
| `stratum`     | factor  | Stratum identifier (27 strata) |
| `psu`         | factor  | PSU identifier (356 PSUs)      |
| `weight`      | numeric | Survey weight                  |
| `has_subsidy` | integer | Binary outcome (0/1)           |

**Usage**

``` r
data(bhf_synthetic_data)
```

**Source**

Generated to mimic the structure of the National Survey of Early Care
and Education (NSECE) listed provider sample.

**Examples**

``` r
data(bhf_synthetic_data)
str(bhf_synthetic_data)
summary(bhf_synthetic_data)

# Domain sample sizes
table(bhf_synthetic_data$state)
```

## Workflow Summary

| Stage | Function | Input | Output |
|----|----|----|----|
| 1\. Compile | [`compile_bhf_model()`](https://joonho112.github.io/bhfvar/reference/compile_bhf_model.md) | — | `stanmodel` |
| 2\. Prepare | [`prepare_bhf_data()`](https://joonho112.github.io/bhfvar/reference/prepare_bhf_data.md) | data.frame | `bhf_data` |
| 3\. Fit | [`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md) | bhf_data, stanmodel | `bhf_fit` |
| 4a. Variance | [`variance_decomposition()`](https://joonho112.github.io/bhfvar/reference/variance_decomposition.md) | bhf_fit | list |
| 4b. Domains | [`domain_estimates()`](https://joonho112.github.io/bhfvar/reference/domain_estimates.md) | bhf_fit | data.frame |
| 4c. Overall | [`overall_estimate()`](https://joonho112.github.io/bhfvar/reference/overall_estimate.md) | bhf_fit | list |
| 4d. LOO | [`log_lik()`](https://joonho112.github.io/bhfvar/reference/log_lik.md) | bhf_fit | matrix |

------------------------------------------------------------------------

*For questions or feedback, please visit the [GitHub
repository](https://github.com/joonho112/bhfvar).*
