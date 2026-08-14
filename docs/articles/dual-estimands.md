# Probability-Scale Estimands A, A\*, and B

## Keep the scales separate

All three named estimands A, A\*, and B are probability-scale
decompositions in the current contract. Latent logit-scale SDs,
variances, and ICCs are returned under `latent`; they are diagnostics,
not aliases for A.

| Quantity | Domain probability | Within component | Intended contrast |
|----|----|----|----|
| A | state-only predictor | Bernoulli | substantive domain distribution |
| A\* | same as A | same as A | A after fixed sampling-noise correction |
| B | full state + stratum + nested-PSU predictor | Bernoulli + design mixture | realized design-composition distribution |

## A

For posterior draw $`d`$,

``` math
p^{A}_{sd}=\operatorname{logit}^{-1}(\alpha_d+u_{sd}).
```

Population shares define the mean, between variance, within variance,
total, and between proportion. The package never silently substitutes
sample shares when exact supplied population shares are requested.

## A\*

A\* subtracts $`\sum_s\pi_s\widehat V_s`$ from A’s between variance and
truncates at zero. It leaves A’s mean and within component unchanged.
Report the correction provenance, `at_boundary`, and `truncated`
summaries. A\* is unavailable when preparation uses
`deattenuation = "none"`.

The correction is conditional on fixed inputs. It does not integrate
over uncertainty in Taylor or externally supplied sampling variances.

## B

B first constructs individual probabilities with all three random-effect
levels, then averages within each domain using within-domain normalized
likelihood weights. Its within-domain variance is

``` math
V_{W,B}=V_{W,\text{binomial}}+V_{W,\text{mixture}}.
```

The mixture term is evaluated by weighted squared deviations, avoiding
the cancellation-prone difference-of-moments form.

## Gaps and API

``` r

vd <- variance_decomposition(fit)
vd$A$summary
vd$A_star$summary
vd$B$summary
vd$gaps$B_minus_A
vd$gaps$A_minus_A_star

domain_estimates(fit, estimand = "A")
domain_estimates(fit, estimand = "B")
```

[`domain_estimates()`](https://joonho112.github.io/bhfvar/reference/domain_estimates.md)
does not expose A*: A* changes a decomposition component, not domain
probabilities. Deprecated `type = "conditional"` maps to A with a
warning. Deprecated `type = "marginal"` errors because it is not a valid
alias for B.

## Reporting rule

Report A\* with its correction source and B with its design-mixture
meaning. Passing oracle tests show that the quantities are computed as
the article defines them; they do not establish the frequentist coverage
of the intervals around them. See
[`vignette("scientific-limitations")`](https://joonho112.github.io/bhfvar/articles/scientific-limitations.md).
