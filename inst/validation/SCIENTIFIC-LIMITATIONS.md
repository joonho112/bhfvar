# Scope and limitations

This file records what the package has been verified to do, what it has not,
and how its output should and should not be interpreted. Read it before using
`bhfvar` for substantive work.

## What has been verified

The implementation is checked against the article's main-text equations and the
Stan program in Appendix B (Lee & Hooper, 2026, *Mathematics* 14(3), 512).

- **Model structure.** Crossed domain and stratum random effects, PSUs nested
  within strata, and the weighted pseudo-likelihood match the article's
  Equations (12) and (13).
- **Centering constraints.** Population-weighted centering of domain effects
  (Equation 16) and within-stratum centering of PSU effects (Equation 17) hold
  to numerical tolerance, and are covered by invariant tests.
- **Estimands.** Probability-scale Estimands A, A*, and B follow Equations
  (19)–(21). A is a design-effect-zero reference standardization. A frozen
  reference implementation, computed independently in base R, reproduces the
  Stan generated quantities.
- **Algebraic identities.** Variance decomposition identities, gap definitions,
  and boundary behaviour are covered by property tests.
- **Package integrity.** The latest full source check recorded 0 errors, 0
  warnings, and 1 `New submission` NOTE; the test suite passed.

These checks establish that the code computes the quantities the article
defines. They do not establish the statistical properties of the resulting
intervals.

## What has not been established

### Interval calibration

**Coverage has been measured in a reduced article-aligned balanced synthetic
outcome-replication design. It is not guaranteed in general or in every tested
condition.**

Across 400 confirmatory fits, pooled central 90% coverage for
`var_between_A`, `var_between_A_star`, `var_between_B`, and
`prop_between_A_star` was 0.900, with a 95% cluster-bootstrap interval of 0.871
to 0.928. The four fixed design cells were 0.930, 0.9175, 0.8925, and 0.860.
The pooled mixture met the preregistered +/-5 percentage-point tolerance around
0.90; this is a tolerance decision, not a claim of exact nominal coverage in
each cell.

The study used a binary outcome, an intercept-only model, 20 domains, 8 strata,
3 PSUs per stratum, 960 observations, and an equal-weight grid crossing two
bundled low/high variance profiles with two weight-informativeness settings. The
article's simulation instead retained a 51-domain empirical survey skeleton and
used a wider scenario grid. The package study is therefore article-aligned but
does not reproduce the article's empirical design, an actual probability
sample, or the restricted-data application.

The fitted target is a survey-weighted pseudo-likelihood, and credible
intervals from an unadjusted pseudo-posterior carry no general frequentist
guarantee. A sandwich-type or design-effect adjustment is not implemented in
this version.

Practical consequence: **treat the intervals as pseudo-posterior credible
intervals.** If your design departs substantially from the conditions above —
in particular if the weights are strongly informative, or if you have far fewer
domains — verify calibration for your own setting, or use replication-based
variance estimation on the point estimates.

### Equivalence and small-effect resolution

The equivalence-curve study is descriptive, not a causal decomposition of what
controls resolution. Its arms were unpaired and changed multiple features. In
particular, the doubled-domain arm also changed observations per cell and the
population-share distribution. The study did not measure a design-effect
covariate, so it does not establish that “design-effect magnitude is dominant”
or that changing sample size or domain count cannot help.

Doubling observations in the selected high-profile comparison reduced the
stored signed relative-gap posterior SD from 0.04480 to 0.03383 (ratio 0.755);
the corresponding folded-SD ratio was 0.822. In the larger curve study, the
median folded-SD ratio for the doubled-observation arm was 0.854. This is real
shrinkage, though less than the simple independent-information benchmark; the
equivalence pass rate did not improve in the unpaired selected cases.

High-profile critical-gap roots required extrapolation outside observed support
and are therefore unidentified. The low-profile root was model-sensitive, so no
universal 1.5% threshold is claimed. Equal-bin curve averages estimate the
designed conditional curve, not the marginal declaration rate under the DGP.

Report the signed gap interval and, if using a binary equivalence rule, its
margin and the observed support. The A-versus-B gap is a
standardization-sensitivity contrast: it measures sensitivity to retaining the
fitted realized design composition and is not specific to informative sampling.

### The article's application

The application reported in the article uses the restricted-use NSECE Level 1
files, which are not distributed with this package. **The package does not
reproduce the article's application estimates**, and the bundled synthetic data
are not a substitute for them.

## Interpretation limits

- **A\* conditions on fixed sampling variances.** Whether supplied directly or
  estimated by Taylor linearization, the $\hat V_s$ are treated as known. Their
  estimation uncertainty is not propagated into the posterior for A*. A large
  zero-boundary mass indicates a boundary-dominated, nonregular diagnostic.
- **A is a reference standardization.** Setting fitted stratum and PSU effects
  to zero does not by itself define a causal or policy intervention.
- **Latent ICC scope.** Latent shares use pre-centering random-effect scale
  parameters, not empirical variances of the realized centered effects.
- **Pseudo log-likelihood is not ordinary log-likelihood.** Aggregating the
  weighted pointwise contributions does not establish ordinary observation- or
  cluster-level LOO/WAIC validity. `log_lik()` exposes both raw and pseudo
  variants so the choice is explicit.
- **Scope.** The model is intercept-only with a binary outcome. Covariates,
  non-binary outcomes, and multi-level extensions are not supported.
- **Scope of the plots.** The plotting helpers (`bhf_plot_variance()`,
  `bhf_plot_icc()`, `bhf_plot_domains()`, `bhf_plot_shrinkage()`,
  `bhf_plot_astar_sensitivity()`) need `ggplot2`, a suggested dependency.
  Without it they fail with a message rather than a broken plot.

## Reporting practice

When reporting results from this package:

1. Report A, A*, and B together, not selectively.
2. State the provenance of $\hat V_s$ (supplied or Taylor) alongside A*.
3. Describe intervals as pseudo-posterior credible intervals.
4. Inspect `fit$diagnostics` before interpreting any fit.

## References

Lee, J., & Hooper, A. (2026). Disentangling signal from noise: A Bayesian
hybrid framework for variance decomposition in complex surveys with post-hoc
domains. *Mathematics*, 14(3), 512. <https://doi.org/10.3390/math14030512>
