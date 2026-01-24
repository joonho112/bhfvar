# Extract Variance Decomposition Results

Extracts and summarizes the variance decomposition from a fitted BHF
model. Provides estimates for multiple estimands: Policy (A),
Descriptive (B), and De-attenuated (A\*).

## Usage

``` r
variance_decomposition(fit, prob = 0.95, print = TRUE)
```

## Arguments

- fit:

  An object of class `bhf_fit` from
  [`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md).

- prob:

  Numeric. Probability for credible intervals. Default is 0.95.

- print:

  Logical. If TRUE, prints a formatted summary. Default is TRUE.

## Value

A list with components:

- logit:

  Variance components on logit scale (Estimand A)

- prob:

  Variance components on probability scale (Estimand B)

- deatten:

  De-attenuated variance components (Estimand A\*)

- summary_table:

  Formatted summary table

## Details

The function extracts three sets of variance components:

- Estimand A (Policy, logit scale):

  Variance decomposition on the latent logit scale. `icc_state` =
  sigma_state^2 / (sigma_state^2 + sigma_psu^2 + pi^2/3)

- Estimand B (Descriptive, probability scale):

  Variance decomposition on the observed probability scale.
  `icc_prob = Var(p_s) / (Var(p_s) + E(p_s(1-p_s)))`

- Estimand A\* (Policy adjusted, de-attenuated):

  Estimand B with finite-sample variance inflation removed.
  `icc_deatten = (Var(p_s) - V_hat) / (Var(p_s) - V_hat + E(p_s(1-p_s)))`

## Interpretation

The ICC (Intraclass Correlation Coefficient) represents the proportion
of total variance attributable to between-domain differences. Higher
values indicate more geographic heterogeneity.

The difference between ICC_B and ICC_A\* represents the amount of
apparent heterogeneity that is actually due to sampling noise rather
than true substantive differences.

## Examples

``` r
if (FALSE) { # \dontrun{
# After fitting
fit <- bhf_fit(prepared_data, model = model)

# Get variance decomposition
vd <- variance_decomposition(fit)

# Access specific components
vd$logit$icc_mean      # ICC on logit scale
vd$prob$icc_mean       # ICC on probability scale
vd$deatten$icc_mean    # De-attenuated ICC
} # }
```
