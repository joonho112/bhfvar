# Fit the Bayesian Hybrid Framework Model

Fits the BHF model using Hamiltonian Monte Carlo via Stan.

## Usage

``` r
bhf_fit(
  data,
  model = NULL,
  chains = 4,
  iter = 2000,
  warmup = floor(iter/2),
  seed = 1234,
  cores = NULL,
  adapt_delta = 0.95,
  max_treedepth = 12,
  refresh = 200,
  ...
)
```

## Arguments

- data:

  An object of class `bhf_data` from
  [`prepare_bhf_data()`](https://joonho112.github.io/bhfvar/reference/prepare_bhf_data.md).

- model:

  A compiled Stan model from
  [`compile_bhf_model()`](https://joonho112.github.io/bhfvar/reference/compile_bhf_model.md).
  If NULL, the function will attempt to compile the model (with a
  warning).

- chains:

  Integer. Number of MCMC chains. Default is 4.

- iter:

  Integer. Total number of iterations per chain (including warmup).
  Default is 2000.

- warmup:

  Integer. Number of warmup iterations per chain. Default is half of
  `iter`.

- seed:

  Integer. Random seed for reproducibility.

- cores:

  Integer. Number of cores for parallel chains. Default is
  `min(chains, parallel::detectCores())`.

- adapt_delta:

  Numeric between 0 and 1. Target acceptance probability. Higher values
  reduce divergent transitions but slow sampling. Default is 0.95.

- max_treedepth:

  Integer. Maximum tree depth for NUTS. Default is 12.

- refresh:

  Integer. How often to print progress (iterations). Set to 0 to
  suppress output.

- ...:

  Additional arguments passed to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).

## Value

An object of class `bhf_fit` containing:

- schema_version:

  Version of the fit-object contract

- contract_id:

  Stable fit-object contract identifier

- package_version:

  Package version that created the fit

- model_sha256:

  SHA-256 fingerprint of the Stan source

- stanfit:

  The stanfit object from rstan

- data:

  The bhf_data object used for fitting

- model:

  The compiled Stan model

- call:

  The function call

- diagnostics:

  MCMC diagnostics summary

- provenance:

  Data, weight, population-share, sampling-variance, prior,
  sampling-control, and session metadata

## Details

The model uses a non-centered parameterization for random effects, which
typically improves sampling efficiency. Key parameters are:

- alpha:

  Global intercept on logit scale

- sigma_state:

  Between-domain standard deviation

- sigma_stratum:

  Between-stratum design standard deviation

- sigma_psu:

  Between-PSU standard deviation (within stratum)

## Convergence Diagnostics

The function automatically checks for:

- Divergent transitions (should be 0)

- Low Rhat values (should be \< 1.05)

- Legacy rstan `n_eff` (should be \> 100; this is not labeled as bulk or
  tail ESS)

- Tree depth saturation

Warnings are issued if any diagnostics suggest convergence problems.

## Typical Runtime

With default settings (4 chains, 2000 iterations):

- Small data (N \< 500): 1-2 minutes

- Medium data (N ~ 3000): 5-10 minutes

- Large data (N \> 10000): 20+ minutes

## Examples

``` r
if (FALSE) { # \dontrun{
# Step 1: Compile model (once per session)
model <- compile_bhf_model()

# Step 2: Prepare data
data(bhf_synthetic_data)
prepared <- prepare_bhf_data(
  bhf_synthetic_data,
  outcome = "has_subsidy",
  domain = "state",
  strata = "stratum",
  psu = "psu",
  weights = "weight"
)

# Step 3: Fit model
fit <- bhf_fit(prepared, model = model, chains = 4, iter = 2000)

# Step 4: Examine results
print(fit)
summary(fit)
} # }
```
