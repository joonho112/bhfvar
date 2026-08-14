# Scope and limitations

This file records what the package has been verified to do, what it has not,
and how its output should and should not be interpreted. Read it before using
`bhfvar` for substantive work.

## What has been verified

The implementation is checked against the model specified in Appendix A of the
accompanying article (Lee & Hooper, 2026, *Mathematics* 14(3), 512).

- **Model structure.** Crossed domain and stratum random effects, PSUs nested
  within strata, and the weighted pseudo-likelihood match the article's
  Equations (12) and (13).
- **Centering constraints.** Population-weighted centering of domain effects
  (Equation 16) and within-stratum centering of PSU effects (Equation 17) hold
  to numerical tolerance, and are covered by invariant tests.
- **Estimands.** Probability-scale Estimands A, A*, and B follow Equations
  (19)–(21). A frozen reference implementation, computed independently in base
  R, reproduces the Stan generated quantities.
- **Algebraic identities.** Variance decomposition identities, gap definitions,
  and boundary behaviour are covered by property tests.
- **Package integrity.** `R CMD check` completes with no errors, warnings, or
  notes, and the test suite passes.

These checks establish that the code computes the quantities the article
defines. They do not establish the statistical properties of the resulting
intervals.

## What has not been established

### Interval calibration

**Frequentist coverage of the posterior intervals has not been established.**

In an internal simulation study, coverage of central 90% intervals was below
nominal for several quantities. This is consistent with theory rather than
surprising: the fitted target is a survey-weighted pseudo-likelihood, and
credible intervals from an unadjusted pseudo-posterior are not guaranteed to
attain nominal frequentist coverage. A sandwich-type or design-effect
adjustment is not implemented in this version.

Practical consequence: **treat the intervals as pseudo-posterior credible
intervals, not as calibrated confidence intervals.** If you need frequentist
guarantees, use replication-based variance estimation on the point estimates,
or wait for a version that implements a variance adjustment.

### Equivalence and small-effect resolution

Assessing whether a domain gap is practically zero requires a posterior that is
narrow relative to the equivalence margin. At the sample sizes used in internal
testing, the posterior for the relative gap was often wider than a 5% margin in
high-variance settings, so such assessments were inconclusive there. Users
should check the width of the posterior against their own margin before
concluding equivalence.

### The article's application

The application reported in the article uses the restricted-use NSECE Level 1
files, which are not distributed with this package. **The package does not
reproduce the article's application estimates**, and the bundled synthetic data
are not a substitute for them.

## Interpretation limits

- **A\* conditions on fixed sampling variances.** Whether supplied directly or
  estimated by Taylor linearization, the $\hat V_s$ are treated as known. Their
  estimation uncertainty is not propagated into the posterior for A*.
- **Pseudo log-likelihood is not ordinary log-likelihood.** Aggregating the
  weighted pointwise contributions does not establish ordinary observation- or
  cluster-level LOO/WAIC validity. `log_lik()` exposes both raw and pseudo
  variants so the choice is explicit.
- **Scope.** The model is intercept-only with a binary outcome. Covariates,
  non-binary outcomes, and multi-level extensions are not supported.
- **No plotting API.** Build graphics from the tidy extractor output.

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
