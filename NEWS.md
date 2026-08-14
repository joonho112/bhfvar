# bhfvar 0.5.0

This version changes the fitted model. Results are not comparable with those
from version 0.3.0.

## Relationship to the accompanying article

The accompanying article (Lee & Hooper, 2026, *Mathematics* 14(3), 512,
<doi:10.3390/math14030512>) cites `bhfvar` version 0.3.0. Version 0.5.0 is the
version that implements the model given in Appendix A of that article.

If you are following the article, use 0.5.0 or later. Version 0.3.0 remains
available under the `v0.3.0` git tag for anyone who needs the exact code that
existed when the article was published.

## Model changes

* Added a stratum random effect. The linear predictor is now
  `alpha + domain + stratum + psu`, matching Equation (13) of the article.
  Previously it was `alpha + domain + psu`.
* Added population-weighted centering of domain effects (Equation 16) and
  within-stratum centering of PSU effects (Equation 17).
* Variance-component priors are now half-t(3, 0, 2.5) as specified in the
  article. A prior-sensitivity selector is available.
* The intercept prior mean and SD are supplied as data, defaulting to the
  design-weighted logit prevalence and SD 0.5.
* Likelihood weights are globally scaled to mean one by default. The legacy
  domain-level scaling remains available as a deprecated, warned path.

## Estimand changes

* Estimands A, A*, and B are computed on the probability scale following
  Equations (19)–(21).
* Estimand B now averages individual probabilities that include the realized
  stratum and PSU effects, and adds the within-domain mixture variance term.
  The previous marginal approximation has been removed.
* A* subtracts the population-weighted mean sampling variance from A's
  between-domain variance and is floored at zero.
* Added signed diagnostic gaps between A, A*, and B.
* Latent-scale SDs and ICCs are reported separately from the probability-scale
  estimands.

## API changes

* `domain_estimates(type = "conditional")` is deprecated in favour of
  `estimand = "A"`. `type = "marginal"` now errors rather than being silently
  treated as B.
* Removed `validate_input_data()` and `validate_stan_data()`. Validation runs
  inside `prepare_bhf_data()` and reports through its result object.
* `prepare_bhf_data()`, `bhf_fit()`, and the extractors return versioned
  objects with explicit provenance fields.
* `log_lik()` exposes raw and pseudo variants explicitly.
* Added de-attenuation modes: Taylor linearization, supplied variances, or
  disabled, with fail-closed validation.

## Documentation

* Added `vignette("migration-0.5.0")` and `vignette("scientific-limitations")`.
* Rewrote all vignettes and the README against the 0.5.0 contracts.
* Added `inst/validation/CLAIM-MATRIX.md` mapping article components to
  implementation and evidence.

## Scope

Posterior intervals are pseudo-posterior credible intervals. Their frequentist
coverage has not been established; see `vignette("scientific-limitations")`.
The article's restricted-data application is not reproduced by this package.

# bhfvar 0.3.0

Initial experimental implementation. This is the version cited in the
accompanying article.
