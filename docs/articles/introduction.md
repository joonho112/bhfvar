# bhfvar: Bayesian Hybrid Framework for Variance Decomposition

## Introduction

The **bhfvar** package implements the Bayesian Hybrid Framework for
variance decomposition in complex surveys with post-hoc
(non-design-based) domains. When researchers use national surveys to
study geographic variation—such as state-level differences in policy
outcomes—they face a fundamental design-analysis mismatch: the survey
was designed for national estimates, not for the sub-national
comparisons being attempted.

This mismatch creates two critical challenges:

1.  **Informative Sampling**: Complex survey designs (stratification,
    clustering) may be correlated with the outcome, confounding
    substantive domain heterogeneity with design artifacts.

2.  **Finite-Sample Variance Inflation**: When domain sample sizes are
    small, sampling noise inflates apparent between-domain variance,
    conflating signal with noise.

The bhfvar package addresses both challenges through a unified Bayesian
framework that combines design-based principles with model-based
inference.

### Installation

``` r
# From GitHub (development version)
devtools::install_github("joonho112/bhfvar")
```

``` r
library(bhfvar)
library(rstan)
```

## The Core Problem: Signal vs. Noise

### Motivating Example: Child Care Subsidy Receipt

Consider the National Survey of Early Care and Education (NSECE), which
uses a complex multi-stage sampling design to achieve national
representativeness. Suppose a researcher wants to estimate the
proportion of home-based child care providers receiving subsidies in
each U.S. state, and quantify how much this proportion varies across
states.

A naive approach would:

1.  Compute the weighted proportion $`\hat{p}_s`$ for each state $`s`$
2.  Calculate the variance of these estimates across states

However, this approach faces two problems:

**Problem 1: Design Artifacts**

``` r
# Illustration of the problem (conceptual)
# State A might have high subsidy rates NOT because of state policy,
# but because it happens to contain PSUs from high-subsidy strata

# The observed variance mixes:
# - True state effects (policy, economy, demographics)
# - Stratum effects (urban/rural, region)
# - PSU effects (local labor markets, provider networks)
```

**Problem 2: Sampling Noise**

``` r
# A state with only 5 sampled providers will have a noisy estimate
# This noise gets counted as "between-state variance"

# E[Var_between_naive] ≈ Var_between_true + Σ π_s × V_s
#                                           ↑
#                           This term inflates the estimate!
```

## The Bayesian Hybrid Framework Solution

### Key Innovations

The bhfvar package implements several key methodological innovations:

1.  **Bayesian Pseudo-Likelihood (BPL)**: Incorporates survey weights
    into the likelihood to achieve design-consistent inference while
    maintaining the benefits of Bayesian multilevel modeling.

2.  **Hybrid GLMM Structure**: Simultaneously models substantive domain
    effects (states) and nuisance design effects (strata, PSUs),
    correctly partitioning variance attributable to each source.

3.  **Dual Estimand Framework**: Distinguishes between:

    - **Estimand A (Policy)**: Substantive variance net of design
      artifacts
    - **Estimand B (Descriptive)**: Total observed variance
    - **Estimand A\*** (De-attenuated): Policy variance corrected for
      finite-sample inflation

4.  **Integrated De-attenuation**: Both implicit (via Bayesian
    shrinkage) and explicit (method-of-moments correction) approaches to
    remove sampling noise from variance estimates.

### The Model

The core model is a weighted generalized linear mixed model:

``` math

\eta_i = \alpha + u_{\text{state}[i]} + u_{\text{psu}[i]}
```

where:

- $`\alpha`$ is the population-weighted global intercept
- $`u_{\text{state}[s]} \sim N(0, \sigma^2_{\text{state}})`$ captures
  substantive state effects
- $`u_{\text{psu}[j]} \sim N(0, \sigma^2_{\text{psu}})`$ captures
  design-induced PSU effects

The outcome $`Y_i \in \{0, 1\}`$ follows:
``` math

Y_i \mid \eta_i \sim \text{Bernoulli}(\text{logit}^{-1}(\eta_i))
```

Survey weights are incorporated via the pseudo-likelihood:
``` math

L_{\text{BPL}} = \prod_{i \in \mathcal{S}} [P(Y_i \mid \eta_i)]^{w_i^*}
```

## Package Capabilities

The bhfvar package provides a complete workflow for variance
decomposition:

### 1. Model Compilation and Data Preparation

``` r
library(bhfvar)

# Step 1: Compile the Stan model (once per session)
model <- compile_bhf_model()

# Step 2: Prepare your data for Stan
data(bhf_synthetic_data)  # Example dataset

prepared <- prepare_bhf_data(
  data = bhf_synthetic_data,
  outcome = "has_subsidy",
  domain = "state",
  strata = "stratum",
  psu = "psu",
  weights = "weight"
)

print(prepared)
```

### 2. Model Fitting

``` r
# Step 3: Fit the Bayesian model
fit <- bhf_fit(
  data = prepared,
  model = model,
  chains = 4,
  iter = 2000
)

print(fit)
```

### 3. Variance Decomposition

``` r
# Step 4: Extract variance decomposition across all three estimands
vd <- variance_decomposition(fit)

# The output shows:
# - Estimand A: ICC on logit scale (policy estimand)
# - Estimand B: ICC on probability scale (descriptive estimand)
# - Estimand A*: De-attenuated ICC (policy adjusted)
```

### 4. Domain Estimates

``` r
# Step 5: Get domain-specific estimates with shrinkage
estimates <- domain_estimates(fit, type = "marginal")

# Each domain gets:
# - Posterior mean and credible interval
# - Reliability (shrinkage factor)
# - Population share
```

## Understanding the Three Estimands

The package computes three distinct variance decomposition estimands:

### Estimand A: Policy (Logit Scale)

``` math

\text{ICC}_A = \frac{\sigma^2_{\text{state}}}{\sigma^2_{\text{state}} + \sigma^2_{\text{psu}} + \pi^2/3}
```

This estimand operates on the latent logit scale and isolates
substantive state variation by factoring out design effects during model
estimation.

### Estimand B: Descriptive (Probability Scale)

``` math

\text{ICC}_B = \frac{\text{Var}(p_s)}{\text{Var}(p_s) + E[p_s(1-p_s)]}
```

This estimand captures total observed variance on the probability scale,
representing what a researcher would see in the data.

### Estimand A\*: Policy Adjusted (De-attenuated)

``` math

\text{ICC}_{A^*} = \frac{\text{Var}(p_s) - \hat{V}}{\text{Var}(p_s) - \hat{V} + E[p_s(1-p_s)]}
```

where $`\hat{V} = \sum_s \pi_s \hat{V}_s`$ is the population-weighted
average sampling variance. This estimand removes finite-sample inflation
from the descriptive variance.

### Interpreting the Differences

| Estimand | Question Answered | When to Use |
|----|----|----|
| A | How much do latent propensities vary across states? | Theoretical modeling |
| B | How much do observed rates vary across states? | Describing the data landscape |
| A\* | How much substantive variation exists after removing noise? | Policy evaluation |

The gap between B and A\* reveals how much of the apparent geographic
variation is actually attributable to sampling noise rather than true
differences.

## Typical Results Pattern

In our motivating application (NSECE data on child care subsidy
receipt), we found:

``` r
# Illustrative results (from simulation studies)
# 
# Estimand A (logit):      ICC = 0.078 [0.032, 0.146]
# Estimand B (prob):       ICC = 0.042 [0.019, 0.076]
# Estimand A* (de-atten):  ICC = 0.006 [0.005, 0.006]
#
# Interpretation:
# - The observed between-state variation (B) is substantial
# - But most of it is sampling noise, not substantive differences
# - After de-attenuation (A*), only ~0.6% of total variance is between states
# - This is a 7-fold reduction from the naive estimate!
```

This pattern—where de-attenuation reveals much smaller substantive
heterogeneity than naively apparent—is common when domains have small
sample sizes and high sampling variance.

## When to Use bhfvar

### Recommended Use Cases

The bhfvar package is particularly valuable for:

- **Sub-national analysis of national surveys**: When using surveys like
  NSECE, NHES, or ECLS to study state-level variation

- **Domain estimation with small samples**: When some domains have few
  observations, making direct estimation unreliable

- **Variance decomposition with complex designs**: When you need to
  separate substantive variation from design artifacts

- **Policy evaluation requiring geographic comparisons**: When the
  question is whether outcomes truly differ across administrative units

### When Other Approaches May Be Appropriate

Consider alternatives when:

- **Domains are design-based**: If the survey was designed with your
  domains as primary sampling units, standard methods may suffice

- **All domains have large samples**: With 100+ observations per domain,
  direct estimation becomes reliable

- **Design is simple random sampling**: Without stratification or
  clustering, the design-based complications disappear

## Road Map to the Vignettes

The bhfvar package includes comprehensive documentation organized into
two tracks:

### Applied Researchers Track

For users who want to apply the package effectively:

| Vignette | Purpose | Reading Time |
|----|----|----|
| [Quick Start](https://joonho112.github.io/bhfvar/articles/quick-start.md) | Your first variance decomposition in 5 minutes | 5 min |
| [Complete Workflow](https://joonho112.github.io/bhfvar/articles/workflow.md) | End-to-end analysis with real data | 30 min |
| [Diagnostics](https://joonho112.github.io/bhfvar/articles/diagnostics.md) | Convergence checking and model validation | 20 min |

### Methodological Researchers Track

For users interested in the mathematical foundations:

| Vignette | Purpose | Reading Time |
|----|----|----|
| [Methodology](https://joonho112.github.io/bhfvar/articles/methodology.md) | The Bayesian Hybrid Framework in detail | 45 min |
| [Dual Estimands](https://joonho112.github.io/bhfvar/articles/dual-estimands.md) | Understanding A, B, and A\* | 30 min |
| [API Reference](https://joonho112.github.io/bhfvar/articles/api-reference.md) | Complete function documentation | Reference |

## Summary

The bhfvar package addresses a fundamental challenge in geographic
policy analysis: how to correctly decompose variance when using complex
surveys for sub-national inference. Key features include:

1.  **Design-consistent inference** through Bayesian Pseudo-Likelihood

2.  **Correct variance partitioning** via the Hybrid GLMM structure

3.  **Multiple estimands** for different research questions (policy vs. 
    descriptive)

4.  **Automatic de-attenuation** to remove finite-sample variance
    inflation

5.  **Domain-specific estimates** with appropriate uncertainty
    quantification

The package is especially valuable in *low-information settings*—where
domain sample sizes are small and sampling variance is high—because it
reveals how much of the apparent geographic variation is true signal
versus noise.

## References

Lee, J., & Hooper, A. (2025). Disentangling signal from noise: A
Bayesian hybrid framework for variance decomposition in complex surveys
with post-hoc domains. *Mathematics* (under review).

Rabe-Hesketh, S., & Skrondal, A. (2006). Multilevel modelling of complex
survey data. *Journal of the Royal Statistical Society: Series A*,
169(4), 805–827.

Savitsky, T. D., & Toth, D. (2016). Bayesian estimation under
informative sampling. *Electronic Journal of Statistics*, 10(1),
1677–1708.

Zeger, S. L., Liang, K.-Y., & Albert, P. S. (1988). Models for
longitudinal data: A generalized estimating equation approach.
*Biometrics*, 44(4), 1049–1060.

------------------------------------------------------------------------

*For questions or feedback about this package, please visit the [GitHub
repository](https://github.com/joonho112/bhfvar).*
