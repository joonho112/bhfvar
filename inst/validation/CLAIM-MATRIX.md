# Article-to-implementation matrix

Each row maps a component of the accompanying article (Lee & Hooper, 2026,
*Mathematics* 14(3), 512) to its implementation and the evidence for it.

`verified` means the implementation was checked against the article's
specification by line-level comparison, a frozen reference oracle, or an
invariant test. It does **not** mean the statistical properties of the
resulting intervals have been established — see `SCIENTIFIC-LIMITATIONS.md`.

## Model

| Component | Article | Status | Implementation | Evidence |
|---|---|---|---|---|
| Weighted pseudo-likelihood | Eq. (12) | verified | Stan `model` block | `test-stan-core.R` |
| Crossed domain/stratum effects | Eq. (13) | verified | Stan `parameters`, `transformed parameters` | `test-stan-core.R` |
| PSUs nested within strata | §3 | verified | `build_nested_psu_index()` | `test-nested-psu-index.R` |
| Non-centered parameterization | §3.3 | verified | three `z` families | `test-stan-core.R` |
| Population-weighted domain centering | Eq. (16) | verified | `transformed parameters` | centering property tests |
| Stratum mean centering | §3 | verified | `transformed parameters` | centering property tests |
| Within-stratum PSU centering | Eq. (17) | verified | `transformed parameters` | centering property tests |
| Variance-component priors half-t(3, 0, 2.5) | §3.4 | verified | Stan `model` block | prior contract tests |
| Intercept prior supplied as data | §3.4 | verified | `prior_alpha_mean`, `prior_alpha_sd` | `test-data-assembly.R` |
| Weight scaling | §3 | verified | `scale_likelihood_weights()`; global mean-one default | `test-weight-contract.R` |

## Estimands

| Component | Article | Status | Implementation | Evidence |
|---|---|---|---|---|
| Estimand A, probability scale | Eq. (20) | verified | `p_state_A`, A decomposition | oracle parity, property tests |
| Estimand A* de-attenuation | Eq. (19) | verified | `var_between_A_star`, floored at zero | de-attenuation tests |
| Estimand B, realized design mixture | Eq. (21) | verified | `p_state_B` via weighted individual probabilities | B oracle parity |
| B binomial and mixture components | Appendix B | verified | `within_binomial_state_B`, `within_mixture_state_B` | generated-quantity tests |
| Signed diagnostic gaps | §4 | verified | `gap_B_minus_A_*`, `gap_A_minus_A_star_*` | property tests |
| Latent-scale SDs and ICCs | §3 | verified | latent generated quantities | `test-stan-core.R` |
| Finite-population domain dispersion | Appendix B | verified | `sd_state_*` quantities | oracle parity |
| Sampling variance $\hat V_s$ | §2.4 | verified | Taylor and supplied adapters | Taylor oracle tests |
| Pointwise log-likelihood | — | verified | `log_lik_raw`, `log_lik_pseudo` | result contract tests |

## Not implemented or out of scope

| Component | Status | Note |
|---|---|---|
| Interval calibration, in general | **not established** | Pooled primary coverage 0.900 (0.871–0.928) in one reduced balanced synthetic design; fixed-cell range 0.860–0.930; see `SCIENTIFIC-LIMITATIONS.md` |
| Sandwich / design-effect variance adjustment | not implemented | Planned |
| $\hat V_s$ uncertainty propagation into A* | not implemented | A* conditions on fixed $\hat V_s$ |
| Ordinary LOO/WAIC validity | not established | Pseudo-likelihood aggregation |
| Restricted-data NSECE application | not reproduced | Requires Level 1 files |
| Covariates, non-binary outcomes | unsupported | Intercept-only binary model |
| Plotting API | supported | `bhf_plot_*()`; needs `ggplot2` (Suggests) |

## Reference

Lee, J., & Hooper, A. (2026). *Mathematics*, 14(3), 512.
<https://doi.org/10.3390/math14030512>
