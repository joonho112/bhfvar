# Methodology and Implementation Contract

## Pseudo-likelihood model

With binary $`y_i`$, the package evaluates

``` math
\sum_i w_i^*\log\Pr(y_i\mid\eta_i), \qquad
\eta_i=\alpha+u_{s[i]}+v_{h[i]}+b_{j[i]}.
```

Domains $`s`$ and strata $`h`$ are crossed. Each flat PSU $`j`$ is
nested in exactly one stratum. The three random-effect SDs are
`sigma_state`, `sigma_stratum`, and `sigma_psu`.

Default likelihood weights are

``` math
w_i^*=w_i\,N/\sum_i w_i,
```

so they remain positive and sum globally to $`N`$.
`weight_scaling = "legacy_d2"` is deprecated and emits a warning; it is
not the default and is not claimed superior or equivalent.

## Centering and priors

- State effects have population-share-weighted mean zero.
- Stratum effects have simple mean zero.
- PSU effects have mean zero within each stratum.
- The baseline state, stratum, and PSU SD priors are positive-half
  Student-t(3, 0, 2.5); approved sensitivity choices change only the
  state SD prior.

These constraints and the generated quantities have passed static and
numeric oracle parity tests. That implementation evidence is distinct
from recovery performance.

## Probability-scale decompositions

For domain probabilities $`p_s`$ and population shares $`\pi_s`$:

``` math
\bar p=\sum_s\pi_s p_s,\quad
V_B=\sum_s\pi_s(p_s-\bar p)^2,\quad
V_W=\sum_s\pi_s p_s(1-p_s),
```

with total $`V_T=V_B+V_W`$ and proportion $`V_B/V_T`$ when defined. A
uses $`p_s=\operatorname{logit}^{-1}(\alpha+u_s)`$. B instead obtains
domain probabilities from observation-level predictors including stratum
and nested PSU effects. Its within component is the sum of a Bernoulli
component and a within-domain mixture component.

A is therefore a design-effect-zero reference standardization, not a
causal or policy intervention. The A-versus-B gap measures sensitivity
to standardizing over the fitted realized design composition; it is not
specific evidence of informative sampling.

A\* keeps A’s mean and within component while applying

``` math
V_{B,A^*}=\max\left(0,V_{B,A}-\sum_s\pi_s\widehat V_s\right).
```

Boundary and truncation flags are part of the output contract.

The latent ICC quantities use the pre-centering random-effect scale
parameters in the model variance formula. They are not empirical
variances of the realized centered effects and do not by themselves
select a random-effect structure.

## Sampling-variance provenance

`deattenuation = "taylor"` estimates domain variances using a one-stage,
stratified PSU design and raw weights. Singleton strata fail explicitly.
`deattenuation = "supplied"` requires an exact-set named nonnegative
vector and a provenance label. `"none"` marks A\* unavailable. In both
computed modes, sampling variances are treated as fixed in the Stan
model; their estimation uncertainty is not propagated.

## Pseudo-posterior limitation

The model is a survey-weighted Bayesian pseudo-likelihood. Posterior
intervals do not automatically inherit ordinary Bayesian coverage.
`log_lik(kind = "pseudo")` exposes weighted contributions, and
PSU/stratum aggregation only sums those contributions; neither operation
establishes standard LOO/WAIC.

## Evidence status

The implementation matches the article’s specification: the model block,
centering constraints, and generated quantities were checked line by
line against the main-text equations and Appendix B Stan program, and a
frozen reference oracle computed independently in base R reproduces the
Stan output.

What this does not establish is frequentist behaviour outside the tested
scope. In a reduced article-aligned balanced synthetic design, pooled
primary 90% coverage was 0.900 (95% cluster interval 0.871 to 0.928),
with fixed-condition rates from 0.860 to 0.930. The pooled result met a
preregistered +/-5 percentage-point tolerance; it was not evidence of
nominal coverage in every condition, of full posterior calibration, or
of the article’s empirical design. A sandwich-type adjustment is not
implemented here. The article’s restricted-data application is not
reproduced.
