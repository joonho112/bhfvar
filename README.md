# bhfvar: Bayesian Hybrid Framework for Variance Decomposition

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

`bhfvar` decomposes variance across post-hoc domains (such as states) in
complex survey data, separating substantive domain variation from artifacts of
the sampling design and from finite-sample noise. It implements the Bayesian
hybrid framework of Lee & Hooper (2026).

## Which version corresponds to the article?

The article cites `bhfvar` **version 0.3.0**. **Version 0.5.0 is the version
that implements the main-text equations and the Stan program in Appendix B of
the article.**

| Article component | 0.3.0 | 0.5.0 |
|---|:-:|:-:|
| Stratum random effect (Eq. 13) | — | yes |
| Population-weighted domain centering (Eq. 16) | — | yes |
| Within-stratum PSU centering (Eq. 17) | — | yes |
| Estimand A on the probability scale (Eq. 20) | — | yes |
| Estimand A* de-attenuation (Eq. 19) | — | yes |
| Estimand B, realized design mixture (Eq. 21) | — | yes |
| Variance-component priors half-t(3, 0, 2.5) | — | yes |

**If you are following the article, use 0.5.0 or later.** Results from 0.3.0
are not comparable. Version 0.3.0 remains available under the `v0.3.0` git tag
for anyone who needs the exact code that existed at publication.

Note that the exported function names are unchanged between the two versions,
so 0.3.0 scripts will still run under 0.5.0 — but they will return different
numbers. See `vignette("migration-0.5.0")`.

## What the model fits

```text
logit(p_i) = alpha + domain[s_i] + stratum[h_i] + psu_within_stratum[j_i]
```

Domain and stratum effects are crossed; PSU identifiers are nested within
strata. Survey weights enter through a pseudo-likelihood and are globally
scaled to mean one by default.

Three estimands are reported on the probability scale:

- **A** — a design-effect-zero reference standardization based on
  `alpha + domain`; it is not by itself a causal or policy intervention;
- **A\*** — an experimental fixed-input diagnostic that subtracts estimated
  domain sampling variance from A's between-domain component;
- **B** — a standardization over each domain's fitted, realized mixture of
  stratum and PSU effects.

The A-versus-B gap is a standardization-sensitivity contrast; it is not specific
to informative sampling. Latent-scale SDs and ICCs are model-scale diagnostics
based on the pre-centering scale parameters and are reported separately.

## Scope

Posterior intervals are **pseudo-posterior credible intervals**. In 400 fits
from a reduced article-aligned balanced synthetic outcome-replication design,
pooled primary 90% coverage was 0.900 with a 95% cluster-bootstrap interval of
0.871 to 0.928. The four fixed-condition rates ranged from 0.860 to 0.930. The
pooled result met the preregistered ±5 percentage-point tolerance; it does not
establish nominal coverage in every condition, for the article's empirical
design, or for other designs. A sandwich-type variance adjustment is not
implemented in this version.

A* conditions on the sampling variances as fixed and does not propagate their
estimation uncertainty. The model is intercept-only with a binary outcome. The
article's restricted-data NSECE application is not reproduced by this package.

See `vignette("scientific-limitations")` before using the package for
substantive work.

## Installation

```r
# install.packages("remotes")
remotes::install_github("joonho112/bhfvar")
```

A working C++ toolchain and `rstan` are required to compile the bundled Stan
model. The model is compiled on demand rather than at install time.

## Workflow

```r
library(bhfvar)

# 1. Prepare the crossed design
data(bhf_synthetic_data)
prepared <- prepare_bhf_data(
  bhf_synthetic_data,
  outcome        = "has_subsidy",
  domain         = "state",
  strata         = "stratum",
  psu            = "psu",
  weights        = "weight",
  weight_scaling = "mean_one",
  deattenuation  = "taylor"
)

# 2. Compile once per session
model <- compile_bhf_model()

# 3. Fit
fit <- bhf_fit(prepared, model = model, chains = 4, iter = 2000, seed = 1234)

# 4. Extract
vd        <- variance_decomposition(fit)
domains_A <- domain_estimates(fit, estimand = "A")
domains_B <- domain_estimates(fit, estimand = "B")
overall_B <- overall_estimate(fit, estimand = "B")
```

Inspect `fit$diagnostics` before interpreting any fit. Report A, A*, and B
together rather than selectively, and state the provenance of the sampling
variances alongside A*.

## Public API

| Stage | Functions |
|---|---|
| Compile | `compile_bhf_model()`, `get_stan_file_path()` |
| Prepare | `prepare_bhf_data()` |
| Fit | `bhf_fit()` |
| Extract | `variance_decomposition()`, `domain_estimates()`, `overall_estimate()`, `log_lik()` |
| Utility | `calc_eff_n()` |

### Plots and diagnostics

| Purpose | Function |
|---|---|
| Probability-scale decomposition across A, A\*, B | `bhf_plot_variance()` |
| Latent scale-parameter shares | `bhf_plot_icc()` |
| Domain estimates with intervals | `bhf_plot_domains()` |
| Model-based versus weighted direct domain estimates | `bhf_plot_shrinkage()` |
| Sensitivity of A\* to the fixed sampling variances | `bhf_astar_sensitivity()`, `bhf_plot_astar_sensitivity()` |

Plotting needs `ggplot2`, a suggested dependency. Every extractor also returns
tidy output if you would rather build the graphics yourself.

## Citation

```r
citation("bhfvar")
```

> Lee, J., & Hooper, A. (2026). Disentangling signal from noise: A Bayesian
> hybrid framework for variance decomposition in complex surveys with post-hoc
> domains. *Mathematics*, 14(3), 512.
> <https://doi.org/10.3390/math14030512>

## License

MIT © JoonHo Lee
