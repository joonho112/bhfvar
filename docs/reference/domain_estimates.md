# Extract Domain-Level Estimates

Extracts posterior summaries of domain-specific probabilities from a
fitted BHF model.

## Usage

``` r
domain_estimates(fit, type = c("marginal", "conditional"), prob = 0.95)
```

## Arguments

- fit:

  An object of class `bhf_fit` from
  [`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md).

- type:

  Character. Type of probability to extract:

  - `"marginal"`: Marginal probabilities (integrating out within-domain
    variation)

  - `"conditional"`: Conditional probabilities (given domain random
    effect)

  Default is "marginal".

- prob:

  Numeric. Probability for credible intervals. Default is 0.95.

## Value

A data frame with columns:

- domain:

  Domain label (from original data)

- domain_id:

  Domain ID (1:S)

- mean:

  Posterior mean

- sd:

  Posterior standard deviation

- q025:

  Lower credible interval bound

- q500:

  Posterior median

- q975:

  Upper credible interval bound

- pop_share:

  Population share of domain

- reliability:

  Reliability/shrinkage factor for domain

## Details

The two types of probabilities differ in their interpretation:

- Marginal probabilities:

  Average probability for a randomly selected individual from the
  domain, integrating over all uncertainty. These are appropriate for
  population-level inference.

- Conditional probabilities:

  Probability given the estimated domain random effect. These represent
  the model's "best guess" for the domain but ignore uncertainty in the
  random effect.

## Examples

``` r
if (FALSE) { # \dontrun{
# After fitting
fit <- bhf_fit(prepared_data, model = model)

# Get domain estimates
estimates <- domain_estimates(fit, type = "marginal")

# View results
head(estimates)

# Plot estimates
library(ggplot2)
ggplot(estimates, aes(x = reorder(domain, mean), y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.2) +
  coord_flip() +
  labs(x = "Domain", y = "Probability")
} # }
```
