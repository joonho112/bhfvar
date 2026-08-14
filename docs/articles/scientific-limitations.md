# Scope and Limitations

This vignette states what the package has been verified to do, what it
has not, and how to interpret its output. Read it before using `bhfvar`
for substantive work.

## What has been verified

The implementation is checked against the model specified in Appendix A
of the accompanying article (Lee & Hooper, 2026).

| Component | Evidence |
|----|----|
| Weighted pseudo-likelihood, crossed domain/stratum effects, nested PSUs | line-level comparison with Appendix A; Stan contract tests |
| Population-weighted domain centering, within-stratum PSU centering | invariant tests to numerical tolerance |
| Estimands A, A\*, B on the probability scale | frozen reference oracle computed independently in base R |
| Variance decomposition identities, gap definitions, boundary behaviour | algebraic property tests |
| Package integrity | `R CMD check` with no errors, warnings, or notes |

These checks establish that the code computes the quantities the article
defines. They say nothing about the statistical properties of the
intervals around those quantities.

## What has not been established

### Interval calibration

**Coverage has been measured for one design family. It is not
established in general.**

A simulation study of 400 fits under the design used in the accompanying
article found central 90% intervals covering their targets at nominal
rates for the primary estimands: pooled coverage 0.900, 95%
cluster-bootstrap interval 0.871 to 0.928. Across all sixteen monitored
quantities it was 0.892 (0.875 to 0.909). That study used a binary
outcome, an intercept-only model, 20 domains, 8 strata, 3 PSUs per
stratum, 960 observations, and weight informativeness of either zero or
moderate.

**This does not establish coverage for other designs.** The fitted
target is a survey-weighted pseudo-likelihood:

``` math
\log L_{\text{BPL}} = \sum_{i \in \mathcal{S}} w_i^* \log P(Y_i \mid \eta_i)
```

Weighting the log-likelihood changes the curvature of the target in a
way that the posterior does not automatically account for, so credible
intervals from an unadjusted pseudo-posterior are not guaranteed to
attain nominal frequentist coverage. A sandwich-type or design-effect
adjustment is **not implemented in this version**.

Practical consequence: treat the intervals as pseudo-posterior credible
intervals. If your design departs substantially from the conditions
above — in particular if the weights are strongly informative, or if you
have far fewer domains — verify calibration for your own setting, or
apply replication-based variance estimation to the point estimates.

### Equivalence and small-effect resolution

Concluding that a domain gap is practically zero requires a posterior
that is narrow relative to your equivalence margin. In internal testing,
the posterior for the relative gap between Estimands A and B was often
wider than a 5% margin when the design effects were large, so such
conclusions were unavailable there. The estimates were not wrong; the
design could not resolve differences that small.

Doubling the number of observations did **not** help. The precision of
the relative gap is governed by the number of domains, not by the number
of observations within them: the gap is a ratio of between-domain
variance components, and those are estimated from the spread across
domains. Adding PSUs or observations sharpens each domain’s mean while
leaving that spread just as hard to separate from noise.

**Check the width of your posterior against your margin before
concluding equivalence, and expect the number of domains — not the
sample size — to be the binding constraint.** A wide posterior inside a
narrow margin is an inconclusive result, not evidence of equivalence.

### The article’s application

The application in the article uses restricted-use NSECE Level 1 files,
which are not distributed with this package. The package **does not
reproduce the article’s application estimates**, and the bundled
synthetic data are not a substitute.

## Interpretation limits

- **A\* conditions on fixed sampling variances.** Whether supplied
  directly or estimated by Taylor linearization, the $`\hat V_s`$ are
  treated as known. Their estimation uncertainty is not propagated into
  the posterior for A\*.
- **Pseudo log-likelihood is not ordinary log-likelihood.** Aggregating
  the weighted pointwise contributions does not establish ordinary
  observation- or cluster-level LOO/WAIC validity.
  [`log_lik()`](https://joonho112.github.io/bhfvar/reference/log_lik.md)
  exposes both variants so the choice is explicit.
- **Scope.** Intercept-only model, binary outcome. Covariates,
  non-binary outcomes, and further levels are not supported.
- **No plotting API.** Build graphics from the tidy extractor output.

## Why clean diagnostics are not enough

Sampler diagnostics answer whether the sampler explored the posterior
properly. Clean diagnostics are a precondition for interpreting a fit,
not evidence about calibration. The two questions are separate: a
sampler can explore a posterior perfectly while that posterior is the
wrong width for frequentist purposes.

## Reporting practice

1.  Report A, A\*, and B together, not selectively.
2.  State the provenance of $`\hat V_s`$ alongside A\*.
3.  Describe intervals as pseudo-posterior credible intervals.
4.  Inspect `fit$diagnostics` before interpreting any fit.

## Reference

Lee, J., & Hooper, A. (2026). Disentangling signal from noise: A
Bayesian hybrid framework for variance decomposition in complex surveys
with post-hoc domains. *Mathematics*, 14(3), 512.
<https://doi.org/10.3390/math14030512>
