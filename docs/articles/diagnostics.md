# Model Diagnostics and Convergence Checking

## Overview

Bayesian inference via MCMC requires verification that the sampler has
converged to the target distribution. This vignette covers:

1.  Essential diagnostic checks
2.  Interpreting warning messages
3.  Common problems and solutions
4.  Advanced diagnostics with bayesplot

## 1. Quick Diagnostic Check

After fitting a model, always check the basic diagnostics:

``` r

library(bhfvar)

# Fit model (assuming model and prepared exist)
fit <- bhf_fit(prepared, model = model, chains = 4, iter = 2000)

# Print includes diagnostics
print(fit)

# Access diagnostics directly
fit$diagnostics$n_divergent    # Should be 0
fit$diagnostics$rhat_max       # Should be < 1.01
fit$diagnostics$ess_min        # Should be > 400
```

## 2. Understanding the Key Diagnostics

### 2.1 Divergent Transitions

**What they are**: Divergent transitions occur when the Hamiltonian
Monte Carlo (HMC) sampler encounters regions of high curvature that it
cannot navigate accurately.

**How to check**:

``` r

cat("Divergent transitions:", fit$diagnostics$n_divergent, "\n")
```

**Interpretation**:

| Count | Interpretation | Action                                   |
|-------|----------------|------------------------------------------|
| 0     | Good           | Proceed with analysis                    |
| 1–10  | Concerning     | Investigate; may need reparameterization |
| \>10  | Problematic    | Do not trust results; fix model          |

**Solutions**:

``` r

# 1. Increase adapt_delta (default is 0.8)
fit <- bhf_fit(prepared, model = model,
               control = list(adapt_delta = 0.95))

# 2. Increase max_treedepth
fit <- bhf_fit(prepared, model = model,
               control = list(adapt_delta = 0.95, max_treedepth = 12))
```

### 2.2 R-hat (Potential Scale Reduction Factor)

**What it is**: R-hat compares variance within chains to variance
between chains. Values near 1 indicate chains have mixed well.

**How to check**:

``` r

cat("Maximum R-hat:", fit$diagnostics$rhat_max, "\n")

# For detailed R-hat by parameter
library(rstan)
summary(fit$stanfit)$summary[, "Rhat"]
```

**Interpretation**:

| R-hat     | Interpretation                             |
|-----------|--------------------------------------------|
| \< 1.01   | Excellent convergence                      |
| 1.01–1.05 | Acceptable, but more iterations might help |
| 1.05–1.1  | Concerning; increase iterations            |
| \> 1.1    | Not converged; do not use results          |

### 2.3 Effective Sample Size (ESS)

**What it is**: ESS estimates the number of independent samples. MCMC
samples are autocorrelated, so ESS \< actual samples.

**How to check**:

``` r

cat("Minimum ESS:", fit$diagnostics$ess_min, "\n")

# Detailed ESS by parameter
summary(fit$stanfit)$summary[, c("n_eff")]
```

**Interpretation**:

| ESS      | Interpretation                           |
|----------|------------------------------------------|
| \> 1000  | Excellent                                |
| 400–1000 | Good for most purposes                   |
| 100–400  | May be adequate for means, not for tails |
| \< 100   | Run longer; do not trust results         |

## 3. Interpreting Warning Messages

### 3.1 Common Warnings and Their Meaning

**“There were X divergent transitions”**

    Warning: There were 5 divergent transitions after warmup.

→ The sampler encountered problematic geometry. Increase `adapt_delta`.

**“R-hat values above 1.05”**

    Warning: Some Rhat values > 1.05, indicating potential convergence issues.

→ Chains haven’t mixed. Run longer or check model specification.

**“Bulk ESS too low”**

    Warning: Bulk Effective Samples Size (ESS) is too low.

→ Posterior means may be unreliable. Increase iterations.

**“Tail ESS too low”**

    Warning: Tail Effective Samples Size (ESS) is too low.

→ Credible interval endpoints may be unreliable. Increase iterations.

### 3.2 Warning Response Guide

``` r

# If you see convergence warnings, try:

# Step 1: More iterations
fit <- bhf_fit(prepared, model = model,
               chains = 4, iter = 4000, warmup = 2000)

# Step 2: If still divergences, increase adapt_delta
fit <- bhf_fit(prepared, model = model,
               chains = 4, iter = 4000,
               control = list(adapt_delta = 0.95))

# Step 3: If still problems, increase max_treedepth
fit <- bhf_fit(prepared, model = model,
               chains = 4, iter = 4000,
               control = list(adapt_delta = 0.99, max_treedepth = 15))
```

## 4. Visual Diagnostics

### 4.1 Trace Plots

Trace plots show parameter values across iterations. Well-mixed chains
look like “fuzzy caterpillars” with no trends.

``` r

library(rstan)

# Trace plots for key parameters
traceplot(fit$stanfit, pars = c("alpha", "sigma_state", "sigma_psu"))
```

**Good trace plot**: Chains overlap, no trends, looks like white noise

**Bad trace plot**: Chains don’t mix, visible trends, stuck periods

### 4.2 Density Plots

``` r

# Density plots comparing chains
stan_dens(fit$stanfit, pars = c("alpha", "sigma_state", "sigma_psu"),
          separate_chains = TRUE)
```

**Good**: Densities from different chains overlap substantially

**Bad**: Densities are very different across chains

### 4.3 Using bayesplot

``` r

library(bayesplot)

# Extract posterior draws
posterior <- as.array(fit$stanfit)

# MCMC diagnostics
mcmc_trace(posterior, pars = c("alpha", "sigma_state", "sigma_psu"))
mcmc_dens_overlay(posterior, pars = c("alpha", "sigma_state", "sigma_psu"))
mcmc_rank_overlay(posterior, pars = c("alpha", "sigma_state", "sigma_psu"))

# Pairs plot to check correlations
mcmc_pairs(posterior, pars = c("alpha", "sigma_state", "sigma_psu"))
```

## 5. Posterior Predictive Checks

### 5.1 Checking Model Fit

The model generates posterior predictive samples (`y_rep`) that can be
compared to the observed data:

``` r

library(bayesplot)

# Extract observed and replicated data
y <- fit$data$stan_data$y
y_rep <- rstan::extract(fit$stanfit, "y_rep")$y_rep

# Density overlay
ppc_dens_overlay(y, y_rep[1:50, ])

# Proportion of 1s
ppc_stat(y, y_rep, stat = "mean")

# Grouped by domain (if desired)
ppc_stat_grouped(y, y_rep, 
                 group = fit$data$stan_data$state_id,
                 stat = "mean")
```

### 5.2 Checking Domain-Level Predictions

``` r

# Compare predicted vs observed domain proportions
estimates <- domain_estimates(fit, type = "marginal")

# Calculate observed proportions
library(dplyr)
observed <- bhf_synthetic_data %>%
  group_by(state) %>%
  summarize(
    obs_prop = weighted.mean(has_subsidy, weight),
    n = n()
  )

# Merge and compare
comparison <- merge(estimates, observed, by.x = "domain", by.y = "state")

plot(comparison$obs_prop, comparison$mean,
     xlab = "Observed Proportion",
     ylab = "Model Estimate",
     main = "Observed vs. Predicted Domain Proportions")
abline(0, 1, col = "red", lty = 2)
```

## 6. Model Comparison with LOO-CV

### 6.1 Leave-One-Out Cross-Validation

``` r

library(loo)

# Extract log-likelihood
ll <- log_lik(fit)

# Compute LOO
loo_result <- loo(ll)
print(loo_result)

# Check Pareto k diagnostics
plot(loo_result)
```

### 6.2 Interpreting LOO Results

**Pareto k values**:

| k       | Interpretation                                 |
|---------|------------------------------------------------|
| \< 0.5  | Good; reliable estimate                        |
| 0.5–0.7 | Marginal; still useful                         |
| 0.7–1.0 | Poor; observation is influential               |
| \> 1.0  | Very bad; consider removing or moment matching |

## 7. Sensitivity Analysis

### 7.1 Prior Sensitivity

Check if results are sensitive to prior choices:

``` r

# Fit with tighter priors on variance components
# (would require modifying the Stan model or using prior arguments)

# Compare key estimates across different prior specifications
# If results change substantially, priors are influential
```

### 7.2 Weight Scaling Sensitivity

``` r

# Compare results with and without weight scaling
# bhfvar uses Method 2 by default; compare with Method 1 or unweighted

# Prepare data without weights (for comparison)
prepared_unweighted <- prepare_bhf_data(
  bhf_synthetic_data, "has_subsidy", "state", "stratum", "psu",
  weights = NULL
)

fit_unweighted <- bhf_fit(prepared_unweighted, model = model)

# Compare ICC estimates
cat("Weighted ICC_B:", variance_decomposition(fit)$prob$icc_mean, "\n")
cat("Unweighted ICC_B:", variance_decomposition(fit_unweighted)$prob$icc_mean, "\n")
```

## 8. Troubleshooting Guide

### 8.1 Problem: Divergent Transitions

**Symptoms**: Warning about divergent transitions; unreliable estimates

**Diagnosis**:

``` r

# Check which iterations diverged
library(rstan)
sampler_params <- get_sampler_params(fit$stanfit, inc_warmup = FALSE)
divergent <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
print(divergent)
```

**Solutions**: 1. Increase `adapt_delta` to 0.95 or 0.99 2.
Reparameterize (bhfvar already uses non-centered) 3. Check for data
issues (extreme weights, sparse domains)

### 8.2 Problem: Poor Mixing

**Symptoms**: R-hat \> 1.05; chains don’t overlap in trace plots

**Diagnosis**:

``` r

# Check which parameters have high R-hat
summ <- summary(fit$stanfit)$summary
high_rhat <- summ[summ[, "Rhat"] > 1.05, , drop = FALSE]
print(rownames(high_rhat))
```

**Solutions**: 1. Run longer (increase `iter`) 2. Increase warmup period
3. Check for multimodality or identification issues

### 8.3 Problem: Low ESS

**Symptoms**: ESS \< 400 for key parameters

**Diagnosis**:

``` r

# Check which parameters have low ESS
summ <- summary(fit$stanfit)$summary
low_ess <- summ[summ[, "n_eff"] < 400, , drop = FALSE]
print(rownames(low_ess))
```

**Solutions**: 1. Run longer 2. Use thinning (set `thin = 2` or higher)
3. Check for slow-mixing parameters

## 9. Diagnostic Checklist

Before reporting results, verify:

**No divergent transitions** (or understood and addressed)

**R-hat \< 1.01** for all parameters of interest

**ESS \> 400** for all parameters of interest

**Trace plots** show good mixing

**Posterior predictive checks** show reasonable fit

**LOO-CV** shows no highly influential observations

**Sensitivity analysis** shows robust results

## 10. Example Complete Diagnostic Workflow

``` r

library(bhfvar)
library(rstan)
library(bayesplot)

# Fit model
model <- compile_bhf_model()
data(bhf_synthetic_data)
prepared <- prepare_bhf_data(
  bhf_synthetic_data, "has_subsidy", "state", "stratum", "psu", "weight"
)

fit <- bhf_fit(prepared, model = model, chains = 4, iter = 2000)

# 1. Check summary diagnostics
print(fit)

# 2. Visual diagnostics
posterior <- as.array(fit$stanfit)

# Trace plots
mcmc_trace(posterior, pars = c("alpha", "sigma_state", "sigma_psu"))

# Density overlays
mcmc_dens_overlay(posterior, pars = c("alpha", "sigma_state", "sigma_psu"))

# 3. Posterior predictive check
y <- prepared$stan_data$y
y_rep <- rstan::extract(fit$stanfit, "y_rep")$y_rep
ppc_dens_overlay(y, y_rep[1:50, ])

# 4. LOO-CV
library(loo)
ll <- log_lik(fit)
loo_result <- loo(ll)
print(loo_result)

# 5. If all checks pass, proceed with inference
variance_decomposition(fit)
```

## Summary

Proper diagnostic checking is essential for reliable Bayesian inference.
The key checks are:

| Diagnostic  | Threshold          | Interpretation                          |
|-------------|--------------------|-----------------------------------------|
| Divergences | 0                  | Sampler navigated geometry successfully |
| R-hat       | \< 1.01            | Chains have mixed                       |
| ESS         | \> 400             | Sufficient independent samples          |
| Trace plots | Fuzzy caterpillars | Visual confirmation of mixing           |
| PPC         | Data in range      | Model captures data features            |
| LOO         | k \< 0.7           | No highly influential observations      |

When diagnostics fail, **do not proceed with inference**. Address the
underlying issues first.

## References

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, 27(5), 1413–1432.

Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., & Gelman, A.
(2019). Visualization in Bayesian workflow. *Journal of the Royal
Statistical Society: Series A*, 182(2), 389–402.

------------------------------------------------------------------------

*For questions or feedback, please visit the [GitHub
repository](https://github.com/joonho112/bhfvar).*
