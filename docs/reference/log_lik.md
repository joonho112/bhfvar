# Extract Log-Likelihood for LOO-CV

Extracts the log-likelihood matrix for use with the loo package for
leave-one-out cross-validation.

## Usage

``` r
log_lik(fit)
```

## Arguments

- fit:

  An object of class `bhf_fit` from
  [`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md).

## Value

A matrix of dimension (n_iterations x N) containing pointwise
log-likelihood values.

## Examples

``` r
if (FALSE) { # \dontrun{
# After fitting
fit <- bhf_fit(prepared_data, model = model)

# Get log-likelihood and compute LOO
library(loo)
ll <- log_lik(fit)
loo_result <- loo(ll)
print(loo_result)
} # }
```
