# Migrating from 0.3.0 to 0.5.0

## Read this first: your code will run, but the numbers will change

Version 0.5.0 keeps the same exported function names as 0.3.0. **A
script written for 0.3.0 will therefore run under 0.5.0 without error —
and return different numbers.**

This is intentional. Version 0.5.0 fits the model given in Appendix A of
the accompanying article; version 0.3.0 fitted a different model. The
difference is not a bug fix in the numerical sense but a change in what
is being estimated.

**Do not compare results across the two versions.** If you have
published or circulated output from 0.3.0, re-run it under 0.5.0 rather
than assuming the values carry over.

## What changed in the model

|  | 0.3.0 | 0.5.0 |
|----|----|----|
| Linear predictor | `alpha + domain + psu` | `alpha + domain + stratum + psu` |
| Stratum random effect | absent | present |
| Domain centering | none | population-weighted (Eq. 16) |
| PSU centering | none | within-stratum (Eq. 17) |
| Variance-component priors | `normal(0, 1)`, `normal(0, 0.5)` | half-t(3, 0, 2.5) |
| Weight scaling default | domain-level | global mean-one |

Adding the stratum effect and the centering constraints changes every
downstream quantity. In particular, between-domain variance in 0.3.0
absorbed stratum-level variation that 0.5.0 attributes to strata.

## What changed in the estimands

|  | 0.3.0 | 0.5.0 |
|----|----|----|
| Estimand A | latent-scale ICC | probability-scale decomposition (Eq. 20) |
| Estimand B | marginal approximation of domain-only probabilities | domain means over realized stratum and PSU effects, plus within-domain mixture variance (Eq. 21) |
| Estimand A\* | derived from the marginal approximation, floored at 0.001 | derived from A’s between-domain variance, floored at 0 (Eq. 19) |

If you previously reported “Estimand B” from 0.3.0, note that the
quantity did not include design effects. The 0.5.0 Estimand B does.

## API changes

### Removed functions

``` r

validate_input_data()   # removed
validate_stan_data()    # removed
```

Validation now runs inside
[`prepare_bhf_data()`](https://joonho112.github.io/bhfvar/reference/prepare_bhf_data.md)
and reports through its result object. If your script called these
directly, delete the calls — the checks still happen.

### Changed arguments

``` r

# 0.3.0
domain_estimates(fit, type = "conditional")
domain_estimates(fit, type = "marginal")

# 0.5.0
domain_estimates(fit, estimand = "A")   # replaces type = "conditional"
domain_estimates(fit, estimand = "B")   # the descriptive estimand
```

`type = "conditional"` is deprecated and warns. **`type = "marginal"`
now errors** rather than being silently treated as B, because the 0.3.0
marginal quantity is not the 0.5.0 Estimand B.

### New arguments

``` r

prepare_bhf_data(
  ...,
  weight_scaling = "mean_one",   # or "legacy_d2" (deprecated, warns)
  deattenuation  = "taylor"      # or "supplied", "none"
)
```

`deattenuation` is fail-closed: supplying an unusable combination errors
rather than silently disabling A\*.

## Updating a 0.3.0 script

``` r

library(bhfvar)

prepared <- prepare_bhf_data(
  my_data,
  outcome        = "y",
  domain         = "state",
  strata         = "stratum",
  psu            = "psu",
  weights        = "weight",
  weight_scaling = "mean_one",   # new
  deattenuation  = "taylor"      # new
)

model <- compile_bhf_model()
fit   <- bhf_fit(prepared, model = model, chains = 4, iter = 2000, seed = 1)

vd <- variance_decomposition(fit)
a  <- domain_estimates(fit, estimand = "A")   # was type = "conditional"
b  <- domain_estimates(fit, estimand = "B")   # was type = "marginal" — different quantity
```

## Reinstalling

If you have 0.3.0 installed, replace it:

``` r

remotes::install_github("joonho112/bhfvar")
packageVersion("bhfvar")   # confirm 0.5.0 or later
```

Check the version explicitly. Because the function names are unchanged,
an un-upgraded installation gives no visible sign that it is fitting the
older model.

## Which version does the article cite?

The accompanying article cites version 0.3.0. Version 0.5.0 is the
version that implements the model in the article’s Appendix A. If you
are following the article, use 0.5.0 or later. Version 0.3.0 remains
available under the `v0.3.0` git tag.

## Scope

Posterior intervals in 0.5.0 are pseudo-posterior credible intervals
whose frequentist coverage has not been established. See
[`vignette("scientific-limitations")`](https://joonho112.github.io/bhfvar/articles/scientific-limitations.md).
