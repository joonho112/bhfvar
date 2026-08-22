# Changelog

## bhfvar 0.6.0

### New features

- Added plotting helpers:
  [`bhf_plot_variance()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_variance.md)
  for the probability-scale decomposition across Estimands A, A\*, and
  B;
  [`bhf_plot_icc()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_icc.md)
  for the latent-scale variance shares;
  [`bhf_plot_domains()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_domains.md)
  for domain estimates with intervals; and
  [`bhf_plot_shrinkage()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_shrinkage.md)
  for model-based versus weighted raw proportions. These cover the
  graphics the article’s recommended workflow asks for. They need
  `ggplot2`, which is a suggested dependency; without it they fail with
  a message pointing to the tidy extractor output.
- Added
  [`bhf_astar_sensitivity()`](https://joonho112.github.io/bhfvar/reference/bhf_astar_sensitivity.md),
  which recomputes Estimand A\* with the supplied sampling variances
  multiplied by a range of factors, and
  [`bhf_plot_astar_sensitivity()`](https://joonho112.github.io/bhfvar/reference/bhf_plot_astar_sensitivity.md)
  to draw the result. A\* conditions on those variances as known
  (article Section 2.4); this shows what that assumption is worth
  without changing it. No refitting is involved.

### Documentation

- Replaced the earlier causal and threshold claims about A-versus-B
  equivalence. The study arms were unpaired and changed multiple design
  features; the high-profile critical gaps required extrapolation
  outside observed support. The public documentation now reports the
  observed curve and width patterns as descriptive, labels high-profile
  critical gaps unidentified, and does not claim that design-effect
  magnitude is the dominant cause or that a universal 1.5% threshold was
  established.
- Clarified that A is a design-effect-zero reference standardization,
  not a causal policy estimand, and that A-versus-B is a
  standardization-sensitivity contrast rather than an
  informative-sampling-specific diagnostic.
- Corrected the article cross-reference: the mathematical proofs are in
  Appendix A; the Stan program is in Appendix B.

### Scope

- No change to the model or to the estimand definitions.

## bhfvar 0.5.1

### Documentation

- Updated the scope and limitations documents to report a simulation
  study of 400 fits. In a reduced article-aligned balanced synthetic
  design, pooled primary 90% coverage was 0.900 with a 95%
  cluster-bootstrap interval of 0.871 to 0.928. Fixed-condition rates
  ranged from 0.860 to 0.930. The pooled result met the preregistered
  +/-5 percentage-point tolerance; it did not establish nominal coverage
  in every condition or reproduce the article’s empirical design.
- The previous text described coverage as below nominal. That
  description came from a smaller study of 24 fits whose cases had been
  selected for extreme configurations, and whose significance test did
  not account for the correlation between quantities computed from the
  same fit. It was not supported.
- The limitations text now states the conditions under which coverage
  was measured, and continues to note that coverage is not established
  for other designs and that no variance adjustment is implemented.
- Earlier 0.5.1 wording named the number of domains as the binding
  equivalence constraint. The 0.6.0 record withdraws that unsupported
  causal attribution.

### Scope

- No change to the model, estimands, or API.

## bhfvar 0.5.0

This version changes the fitted model. Results are not comparable with
those from version 0.3.0.

### Relationship to the accompanying article

The accompanying article (Lee & Hooper, 2026, *Mathematics* 14(3), 512,
<doi:10.3390/math14030512>) cites `bhfvar` version 0.3.0. Version 0.5.0
is the version that implements the main-text equations and Stan program
in Appendix B of that article.

If you are following the article, use 0.5.0 or later. Version 0.3.0
remains available under the `v0.3.0` git tag for anyone who needs the
exact code that existed when the article was published.

### Model changes

- Added a stratum random effect. The linear predictor is now
  `alpha + domain + stratum + psu`, matching Equation (13) of the
  article. Previously it was `alpha + domain + psu`.
- Added population-weighted centering of domain effects (Equation 16)
  and within-stratum centering of PSU effects (Equation 17).
- Variance-component priors are now half-t(3, 0, 2.5) as specified in
  the article. A prior-sensitivity selector is available.
- The intercept prior mean and SD are supplied as data, defaulting to
  the design-weighted logit prevalence and SD 0.5.
- Likelihood weights are globally scaled to mean one by default. The
  legacy domain-level scaling remains available as a deprecated, warned
  path.

### Estimand changes

- Estimands A, A\*, and B are computed on the probability scale
  following Equations (19)–(21).
- Estimand B now averages individual probabilities that include the
  realized stratum and PSU effects, and adds the within-domain mixture
  variance term. The previous marginal approximation has been removed.
- A\* subtracts the population-weighted mean sampling variance from A’s
  between-domain variance and is floored at zero.
- Added signed diagnostic gaps between A, A\*, and B.
- Latent-scale SDs and ICCs are reported separately from the
  probability-scale estimands.

### API changes

- `domain_estimates(type = "conditional")` is deprecated in favour of
  `estimand = "A"`. `type = "marginal"` now errors rather than being
  silently treated as B.
- Removed
  [`validate_input_data()`](https://joonho112.github.io/bhfvar/reference/validate_input_data.md)
  and
  [`validate_stan_data()`](https://joonho112.github.io/bhfvar/reference/validate_stan_data.md).
  Validation runs inside
  [`prepare_bhf_data()`](https://joonho112.github.io/bhfvar/reference/prepare_bhf_data.md)
  and reports through its result object.
- [`prepare_bhf_data()`](https://joonho112.github.io/bhfvar/reference/prepare_bhf_data.md),
  [`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md),
  and the extractors return versioned objects with explicit provenance
  fields.
- [`log_lik()`](https://joonho112.github.io/bhfvar/reference/log_lik.md)
  exposes raw and pseudo variants explicitly.
- Added de-attenuation modes: Taylor linearization, supplied variances,
  or disabled, with fail-closed validation.

### Documentation

- Added
  [`vignette("migration-0.5.0")`](https://joonho112.github.io/bhfvar/articles/migration-0.5.0.md)
  and
  [`vignette("scientific-limitations")`](https://joonho112.github.io/bhfvar/articles/scientific-limitations.md).
- Rewrote all vignettes and the README against the 0.5.0 contracts.
- Added `inst/validation/CLAIM-MATRIX.md` mapping article components to
  implementation and evidence.

### Scope

Posterior intervals are pseudo-posterior credible intervals. Their
coverage has been measured only in a reduced balanced synthetic design
and is not guaranteed more generally; see
[`vignette("scientific-limitations")`](https://joonho112.github.io/bhfvar/articles/scientific-limitations.md).
The article’s restricted-data application is not reproduced by this
package.

## bhfvar 0.3.0

Initial experimental implementation. This is the version cited in the
accompanying article.
