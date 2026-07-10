# bhfvar: Bayesian Hybrid Framework for Variance Decomposition

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Overview

**bhfvar** (Bayesian Hybrid Framework for VARiance decomposition) implements the methodology from Lee & Hooper (2026) for variance decomposition in complex surveys with post-hoc domains.

The package provides tools for:

- **Separating signal from noise**: Distinguish true geographic variation from sampling artifacts
- **Dual estimand framework**: Policy-relevant (latent) vs. Descriptive (observed) variance decomposition
- **De-attenuation**: Correct for finite-sample variance inflation

## Installation

```r
# Install from GitHub (requires devtools)
# install.packages("devtools")
devtools::install_github("joonho112/bhfvar")
```

### Prerequisites

This package requires a working C++ toolchain for Stan compilation:

- **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
- **macOS**: Run `xcode-select --install` in Terminal
- **Linux**: Install build-essential (e.g., `sudo apt install build-essential`)

## Quick Start

```r
library(bhfvar)

# Step 1: Compile the Stan model (once per session)
model <- compile_bhf_model()

# Step 2: Load example data
data(bhf_synthetic_data)

# Step 3: Prepare data for Stan
prepared <- prepare_bhf_data(
  bhf_synthetic_data,
  outcome = "has_subsidy",
  domain = "state",
  strata = "stratum",
  psu = "psu",
  weights = "weight"
)

# Step 4: Fit the model
fit <- bhf_fit(prepared, model = model, chains = 4, iter = 2000)

# Step 5: Extract results
variance_decomposition(fit)
domain_estimates(fit, type = "marginal")
```

## Key Features

### Three Estimands

| Estimand | Scale | Interpretation |
|----------|-------|----------------|
| **A** (Policy) | Logit | Latent heterogeneity in underlying propensity |
| **B** (Descriptive) | Probability | Observed heterogeneity in rates |
| **A*** (Adjusted) | Probability | True heterogeneity after removing sampling noise |

### Defensive Programming Approach

Unlike packages that pre-compile Stan models during installation, `bhfvar` uses a manual compilation approach:

```r
# Compile once per R session
model <- compile_bhf_model()
```

**Benefits:**
- More stable across R/Stan/C++ version combinations
- Clear error messages when compilation fails
- Full visibility into the compilation process
- No hidden cached objects causing mysterious errors

### Comprehensive Diagnostics

```r
# Model diagnostics
fit$diagnostics

# Convergence checks
print(fit)

# Posterior predictive checks (requires bayesplot)
library(bayesplot)
y_rep <- rstan::extract(fit$stanfit, "y_rep")$y_rep
bayesplot::ppc_dens_overlay(bhf_synthetic_data$has_subsidy, y_rep[1:50, ])
```

## Methodology

The package implements a Bayesian multilevel model with:

1. **Pseudo-likelihood weighting**: Incorporates survey weights via Method D2 scaling
2. **Non-centered parameterization**: Improves MCMC sampling efficiency
3. **Marginal scaling**: Converts logit-scale effects to probability scale using Zeger et al. (1988) approximation
4. **De-attenuation**: Removes finite-sample variance inflation from ICC estimates

### Model Structure

```
y_i ~ Bernoulli(p_i)
logit(p_i) = alpha + u_state[s_i] + u_psu[j_i]

u_state ~ Normal(0, sigma_state)
u_psu ~ Normal(0, sigma_psu)
```

### Variance Decomposition

**Probability scale (Estimand B):**
```
Var_between = Sum_s{ pi_s * (p_s - p_bar)^2 }
Var_within = Sum_s{ pi_s * p_s * (1 - p_s) }
ICC_B = Var_between / (Var_between + Var_within)
```

**De-attenuated (Estimand A*):**
```
Var_between_adj = max(0, Var_between - V_hat)
ICC_A* = Var_between_adj / (Var_between_adj + Var_within)
```

where `V_hat` is the estimated sampling variance due to survey design.

## Development Scripts

The `dev/` folder contains step-by-step scripts for understanding and testing the package:

| Script | Purpose |
|--------|---------|
| `run_workflow.R` | **Master script** - runs complete workflow from source |
| `00_verify_environment.R` | Check R and Stan setup |
| `01_compile_model.R` | Manual Stan model compilation |
| `02_prepare_data.R` | Data preparation and validation |
| `03_fit_model.R` | Model fitting with diagnostics |
| `04_extract_results.R` | Result extraction and visualization |

### Quick Test (from source)

```r
# From package root directory:
setwd("path/to/bhfvar")
source("dev/run_workflow.R")
```

## Model Structure (2-Level)

This package implements a **2-level hierarchical model**:

```
Level 2: States (S domains)
  - u_state[s] ~ Normal(0, sigma_state)
  
Level 1: PSUs (J clusters nested in strata)  
  - u_psu[j] ~ Normal(0, sigma_psu)
```

**Key parameters:**
- `alpha`: Global intercept (logit scale)
- `sigma_state`: Between-state standard deviation
- `sigma_psu`: Between-PSU standard deviation

**Generated quantities:**
- `icc_state`: ICC on logit scale (Estimand A)
- `icc_prob`: ICC on probability scale (Estimand B)
- `icc_deatten`: De-attenuated ICC (Estimand A*)

## Citation

If you use this package, please cite the companion paper:

```
Lee, J., & Hooper, A. (2026). Disentangling signal from noise: A Bayesian
hybrid framework for variance decomposition in complex surveys with
post-hoc domains. Mathematics, 14(3), 512.
https://doi.org/10.3390/math14030512
```

Run `citation("bhfvar")` in R for BibTeX entries for both the paper and the
package itself.

## License
MIT © JoonHo Lee

## Contact

- **Author**: JoonHo Lee
- **Email**: jlee296@ua.edu
- **Institution**: The University of Alabama
