# Understanding the Dual Estimand Framework: A, B, and A\*

## Overview

When analyzing geographic variation in complex surveys, the question
“How much do outcomes vary across domains?” has multiple valid
interpretations. The bhfvar package implements the **Dual Estimand
Framework** that distinguishes between three conceptually distinct
quantities:

| Estimand | Name            | Scale       | Purpose                            |
|----------|-----------------|-------------|------------------------------------|
| A        | Policy          | Logit       | Latent propensity variation        |
| B        | Descriptive     | Probability | Observed rate variation            |
| A\*      | Policy Adjusted | Probability | Substantive variation net of noise |

This vignette explains each estimand, when to use it, and how to
interpret the differences between them.

## 1. The Fundamental Ambiguity

### 1.1 What Does “Between-State Variance” Mean?

Consider the fitted Hybrid GLMM:
``` math
\eta_i = \alpha + u_{\text{state}[i]} + u_{\text{psu}[i]}
```

After estimation, we have variance components
$`\sigma^2_{\text{state}}`$ and $`\sigma^2_{\text{psu}}`$. But to
compute a population-level ICC on the probability scale, we must decide
how to handle the PSU effects.

This leads to two interpretations:

**Substantive Interpretation (Policy Question)** \> How much would
outcomes vary across states if we could conceptually \> neutralize the
sampling design?

**Descriptive Interpretation (Data Landscape)** \> How much do outcomes
actually vary across states in the population, \> given their specific
composition of PSUs?

These are different questions requiring different estimands.

### 1.2 The Role of Design Effects

If the sampling design is **informative**—PSUs correlate with the
outcome—then PSU effects are distributed unevenly across states. A state
might appear to have high subsidy rates not because of state policy, but
because it contains more high-subsidy PSUs.

The Hybrid GLMM separates these effects during estimation. But when
*summarizing* results on the probability scale, we must choose whether
to include or exclude them.

## 2. Estimand A: The Policy Estimand (Logit Scale)

### 2.1 Definition

Estimand A operates on the latent logit scale:
``` math
\text{ICC}_A = \frac{\sigma^2_{\text{state}}}{\sigma^2_{\text{state}} + 
\sigma^2_{\text{psu}} + \pi^2/3}
```

where $`\pi^2/3 \approx 3.29`$ is the level-1 variance of the standard
logistic distribution.

### 2.2 Interpretation

Estimand A answers: *What proportion of latent propensity variation is
between states?*

- The numerator $`\sigma^2_{\text{state}}`$ is substantive state
  variation
- The denominator includes all sources: states, PSUs, and individual
  variation
- The ICC measures how much knowing someone’s state tells you about
  their latent propensity

### 2.3 When to Use

Use Estimand A when:

- You want a summary on the modeling scale (logit)
- You’re comparing with other random effects logistic models
- You need a quick assessment of relative variance magnitudes

### 2.4 Accessing in bhfvar

``` r

vd <- variance_decomposition(fit)

# Logit-scale results
cat("Between-state variance (logit):", vd$logit$var_between_mean, "\n")
cat("Within-state variance (logit):", vd$logit$var_within_mean, "\n")
cat("ICC (Estimand A):", vd$logit$icc_mean, "\n")
cat("95% CI:", vd$logit$icc_q025, "-", vd$logit$icc_q975, "\n")
```

## 3. Estimand B: The Descriptive Estimand (Probability Scale)

### 3.1 Definition

Estimand B operates on the probability scale:
``` math
\text{ICC}_B = \frac{\sum_{s} \pi_s (p_s - \bar{p})^2}{\sum_{s} \pi_s (p_s - \bar{p})^2 + \sum_{s} \pi_s \cdot p_s(1-p_s)}
```

where:

- $`p_s = \text{logit}^{-1}(c \cdot (\alpha + u_{\text{state}[s]}))`$ is
  the marginal probability for state $`s`$
- $`c`$ is the Zeger adjustment factor for integrating out PSU variation
- $`\bar{p} = \sum_s \pi_s p_s`$ is the population mean probability

### 3.2 Interpretation

Estimand B answers: *What proportion of total variance in rates is
between states?*

- The numerator is the weighted variance of state probabilities
- The denominator adds expected within-state (Bernoulli) variance
- This represents what you would observe in the actual population

### 3.3 Key Properties

1.  **Observable scale**: Directly interpretable as variation in rates
2.  **Includes all sources**: Design effects influence $`p_s`$ through
    the marginal integration
3.  **Subject to sampling noise**: Small-sample estimates of $`p_s`$
    inflate the numerator

### 3.4 When to Use

Use Estimand B when:

- You want to describe the actual landscape of rates across states
- You’re communicating with non-technical audiences
- The research question is purely descriptive (not causal/policy)

### 3.5 Accessing in bhfvar

``` r

vd <- variance_decomposition(fit)

# Probability-scale results
cat("Between-state variance (prob):", vd$prob$var_between_mean, "\n")
cat("Within-state variance (prob):", vd$prob$var_within_mean, "\n")
cat("ICC (Estimand B):", vd$prob$icc_mean, "\n")
cat("95% CI:", vd$prob$icc_q025, "-", vd$prob$icc_q975, "\n")
```

## 4. Estimand A\*: The Policy Adjusted Estimand (De-attenuated)

### 4.1 The Problem with Estimand B

While Estimand B is interpretable, it conflates two sources of
variation:

1.  **True substantive differences** across states
2.  **Sampling noise** from estimating state proportions

When domain sample sizes are small, sampling variance can be
substantial:
``` math
\hat{p}_s = p_s + e_s, \quad \text{Var}(e_s) = V_s
```

The naive between-state variance estimate is inflated:
``` math
E[\text{Var}(\hat{p}_s)] \approx \text{Var}(p_s) + \sum_s \pi_s V_s
```

### 4.2 Definition

Estimand A\* corrects for this inflation:
``` math
\text{Var}_{\text{between}}^{*} = \max\left(0, \text{Var}_{\text{between}} - 
\sum_{s=1}^S \pi_s \hat{V}_s\right)
```

where $`\hat{V}_s`$ is the design-based sampling variance for state
$`s`$.

The de-attenuated ICC is:
``` math
\text{ICC}_{A^*} = \frac{\text{Var}_{\text{between}}^{*}}{\text{Var}_{\text{between}}^{*} + \text{Var}_{\text{within}}}
```

### 4.3 Interpretation

Estimand A\* answers: *What proportion of substantive variance is
between states, after removing sampling noise?*

- It targets the same quantity as B but corrects for finite-sample bias
- When domains are well-sampled, A\* ≈ B
- When domains have small samples, A\* \< B (potentially much smaller)

### 4.4 When to Use

Use Estimand A\* when:

- You want to know how much *true* heterogeneity exists
- Domain sample sizes are uneven or small
- You’re making policy-relevant comparisons
- You need to separate signal from noise

### 4.5 Accessing in bhfvar

``` r

vd <- variance_decomposition(fit)

# De-attenuated results
cat("Between-state variance (de-atten):", vd$deatten$var_between_mean, "\n")
cat("ICC (Estimand A*):", vd$deatten$icc_mean, "\n")
cat("95% CI:", vd$deatten$icc_q025, "-", vd$deatten$icc_q975, "\n")

# Compare B vs A* to see how much is noise
cat("\nNoise assessment:\n")
cat("ICC_B / ICC_A*:", round(vd$prob$icc_mean / vd$deatten$icc_mean, 2), "\n")
```

## 5. Comparing the Estimands

### 5.1 Typical Relationships

In most applications, you’ll find:
``` math
\text{ICC}_{A^*} < \text{ICC}_B \leq \text{ICC}_A
```

The gap between B and A\* reveals how much apparent heterogeneity is
actually sampling noise.

### 5.2 Worked Example

``` r

library(bhfvar)
library(rstan)

# Fit model
model <- compile_bhf_model()
data(bhf_synthetic_data)
prepared <- prepare_bhf_data(
  bhf_synthetic_data, "has_subsidy", "state", "stratum", "psu", "weight"
)
fit <- bhf_fit(prepared, model = model, chains = 4, iter = 2000)

# Compare all three estimands
vd <- variance_decomposition(fit)

comparison <- data.frame(
  Estimand = c("A (Logit)", "B (Probability)", "A* (De-attenuated)"),
  ICC = c(vd$logit$icc_mean, vd$prob$icc_mean, vd$deatten$icc_mean),
  Lower = c(vd$logit$icc_q025, vd$prob$icc_q025, vd$deatten$icc_q025),
  Upper = c(vd$logit$icc_q975, vd$prob$icc_q975, vd$deatten$icc_q975)
)

print(comparison)
```

### 5.3 Interpretation Guide

| Scenario   | Interpretation                                     |
|------------|----------------------------------------------------|
| B ≈ A\*    | Sample sizes are adequate; little noise inflation  |
| B \>\> A\* | Much of apparent variation is noise; small samples |
| A \> B     | Expected due to logit vs. probability scale        |
| A\* ≈ 0    | Very little true between-state variation           |

## 6. Visualizing the Estimands

### 6.1 ICC Comparison Plot

``` r

library(ggplot2)

# Create comparison data frame
icc_data <- data.frame(
  Estimand = factor(c("A (Logit)", "B (Prob)", "A* (De-atten)"),
                    levels = c("A (Logit)", "B (Prob)", "A* (De-atten)")),
  ICC = c(vd$logit$icc_mean, vd$prob$icc_mean, vd$deatten$icc_mean),
  Lower = c(vd$logit$icc_q025, vd$prob$icc_q025, vd$deatten$icc_q025),
  Upper = c(vd$logit$icc_q975, vd$prob$icc_q975, vd$deatten$icc_q975)
)

ggplot(icc_data, aes(x = Estimand, y = ICC)) +
  geom_point(size = 4, color = "#377EB8") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "#377EB8") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "ICC Comparison Across Estimands",
    subtitle = "With 95% credible intervals",
    y = "Intraclass Correlation Coefficient",
    x = NULL
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11))
```

### 6.2 Domain-Level View

``` r

# Get domain estimates
estimates <- domain_estimates(fit, type = "marginal")

# Sort by estimate
estimates <- estimates[order(estimates$mean), ]
estimates$rank <- 1:nrow(estimates)

ggplot(estimates, aes(x = rank, y = mean)) +
  geom_point(aes(size = pop_share), color = "#377EB8", alpha = 0.7) +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0, alpha = 0.5) +
  geom_hline(aes(yintercept = overall_estimate(fit)$mean), 
             linetype = "dashed", color = "red") +
  labs(
    title = "Domain Estimates (Caterpillar Plot)",
    subtitle = "Sorted by posterior mean; point size = population share",
    x = "Domain Rank",
    y = "Probability",
    size = "Pop Share"
  ) +
  theme_minimal()
```

## 7. Decision Framework

### 7.1 Which Estimand for Which Question?

| Research Question                          | Recommended Estimand |
|--------------------------------------------|----------------------|
| “How much do latent propensities vary?”    | A                    |
| “What does the data landscape look like?”  | B                    |
| “How much true heterogeneity exists?”      | A\*                  |
| “Are state differences meaningful?”        | Compare B vs. A\*    |
| “Should we target interventions by state?” | A\*                  |

### 7.2 Reporting Recommendations

**For policy research**: Report A\* as the primary estimand, with B as
context for what the “naive” analysis would show.

**For descriptive research**: Report B, but note that it may be inflated
by sampling noise if domain samples are small.

**For methods papers**: Report all three with discussion of their
relationships.

### 7.3 Example Write-Up

> “The observed between-state variance in subsidy receipt (ICC$`_B`$ =
> 0.042, 95% CI: 0.019–0.076) substantially exceeds the de-attenuated
> estimate (ICC$`_{A^*}`$ = 0.006, 95% CI: 0.005–0.006). This seven-fold
> difference indicates that most of the apparent geographic variation is
> attributable to sampling noise rather than true substantive
> differences. After de-attenuation, only approximately 0.6% of total
> variance in subsidy receipt lies between states.”

## 8. Technical Details

### 8.1 De-attenuation Within the Posterior

The de-attenuation is performed **at each MCMC iteration**:

1.  Compute $`\text{Var}_{\text{between}}^{(t)}`$ for iteration $`t`$
2.  Subtract $`\sum_s \pi_s \hat{V}_s`$ (fixed, pre-computed from data)
3.  Apply floor at small positive value (0.001) to avoid numerical
    issues
4.  Compute $`\text{ICC}_{A^*}^{(t)}`$

This propagates uncertainty from the model parameters into the
de-attenuated estimates.

### 8.2 Sampling Variance Estimation

The $`\hat{V}_s`$ are computed via Taylor linearization in
[`prepare_bhf_data()`](https://joonho112.github.io/bhfvar/reference/prepare_bhf_data.md):

``` r

# The prepared data includes sampling variances
prepared$sampling_variances  # Vector of V_s

# Population-weighted average
vhat_mean <- sum(prepared$stan_data$w_state_pop_share * 
                 prepared$stan_data$vhat_state)
```

### 8.3 Edge Cases

**All domains well-sampled**: A\* ≈ B; de-attenuation has little effect

**Some domains very sparse**: A\* may be much smaller than B; watch for
domains where $`n_s < 10`$

**Negative de-attenuated variance**: Floored at 0.001; indicates true
between-domain variance is negligible

## Summary

The Dual Estimand Framework provides a principled approach to variance
decomposition that acknowledges the distinct questions researchers may
ask:

| Estimand | Scale       | Includes Noise?     | Use Case      |
|----------|-------------|---------------------|---------------|
| A        | Logit       | Corrected by design | Theoretical   |
| B        | Probability | Yes                 | Descriptive   |
| A\*      | Probability | Corrected           | Policy/Causal |

The key insight is that **apparent heterogeneity ≠ true heterogeneity**
when domain samples are small. By comparing B and A\*, researchers can
quantify how much of the observed variation is signal versus noise.

## References

Lee, J., & Hooper, A. (2026). Disentangling signal from noise: A
Bayesian hybrid framework for variance decomposition in complex surveys
with post-hoc domains. *Mathematics*, 14(3), 512.
<https://doi.org/10.3390/math14030512>

Zeger, S. L., Liang, K.-Y., & Albert, P. S. (1988). Models for
longitudinal data: A generalized estimating equation approach.
*Biometrics*, 44(4), 1049–1060.

------------------------------------------------------------------------

*For questions or feedback, please visit the [GitHub
repository](https://github.com/joonho112/bhfvar).*
