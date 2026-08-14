# Synthetic Survey Data for BHF Package Examples

A legacy synthetic dataset retained for examples and 0.3.0 migration. It
mimics a complex survey structure but contains no restricted NSECE data.
It is not the scientific-validation fixture and must not be used to
claim reproduction of the article's restricted-data results.

## Usage

``` r
bhf_synthetic_data
```

## Format

A data frame with 1,598 observations and 5 variables:

- state:

  State identifier (character). 50 unique states.

- stratum:

  Stratum identifier (character). 27 unique strata.

- psu:

  Primary Sampling Unit identifier (character). Nested within strata.

- weight:

  Survey sampling weight (numeric). Reflects unequal probability
  sampling.

- has_subsidy:

  Binary outcome variable (integer, 0/1). Indicates subsidy receipt.

## Source

Synthetic data generated for the bhfvar package. Not derived from actual
survey responses.

## Details

The data were generated with the following characteristics:

- Overall proportion of has_subsidy approximately 0.27

- State-level random effects with SD approximately 0.5 (logit scale)

- PSU-level random effects with SD approximately 0.3 (logit scale)

- Realistic weight variation reflecting complex survey design

- State population shares varying from small to large states

This synthetic data is suitable for:

- Testing package functions

- Running package examples

- Learning the BHF methodology

- Verifying installation

## Provenance

The exact historical generator is not part of the current source
contract. The package preserves this data object for compatibility;
scientific validation uses separately versioned fixtures and frozen
provenance.

## Examples

``` r
# Load the data
data(bhf_synthetic_data)

# View structure
str(bhf_synthetic_data)
#> 'data.frame':    1598 obs. of  5 variables:
#>  $ state      : chr  "State_01" "State_01" "State_01" "State_01" ...
#>  $ stratum    : chr  "Stratum_02" "Stratum_02" "Stratum_02" "Stratum_02" ...
#>  $ psu        : chr  "PSU_033" "PSU_042" "PSU_022" "PSU_022" ...
#>  $ weight     : num  73 62.1 63.8 145.4 241.7 ...
#>  $ has_subsidy: int  1 0 1 0 1 0 0 1 1 0 ...

# Summary statistics
table(bhf_synthetic_data$has_subsidy)
#> 
#>    0    1 
#> 1142  456 
length(unique(bhf_synthetic_data$state))
#> [1] 50

# State-level proportions
aggregate(has_subsidy ~ state, data = bhf_synthetic_data, FUN = mean)
#>       state has_subsidy
#> 1  State_01  0.38461538
#> 2  State_02  0.18181818
#> 3  State_03  0.23333333
#> 4  State_04  0.46202532
#> 5  State_05  0.35714286
#> 6  State_06  0.33333333
#> 7  State_07  0.60000000
#> 8  State_08  0.53846154
#> 9  State_09  0.12000000
#> 10 State_10  0.34000000
#> 11 State_11  0.19607843
#> 12 State_12  0.50000000
#> 13 State_13  0.31578947
#> 14 State_14  0.06250000
#> 15 State_15  0.15217391
#> 16 State_16  0.18750000
#> 17 State_17  0.10344828
#> 18 State_18  0.11764706
#> 19 State_19  0.26923077
#> 20 State_20  0.35000000
#> 21 State_21  0.36792453
#> 22 State_22  0.15000000
#> 23 State_23  0.31578947
#> 24 State_24  0.18181818
#> 25 State_25  0.00000000
#> 26 State_26  0.43181818
#> 27 State_27  0.05000000
#> 28 State_28  0.37500000
#> 29 State_29  0.20000000
#> 30 State_30  0.00000000
#> 31 State_31  0.30769231
#> 32 State_32  0.50000000
#> 33 State_33  0.27272727
#> 34 State_34  0.00000000
#> 35 State_35  0.40000000
#> 36 State_36  0.17241379
#> 37 State_37  0.25000000
#> 38 State_38  0.07142857
#> 39 State_39  0.14583333
#> 40 State_40  0.33928571
#> 41 State_41  0.43478261
#> 42 State_42  0.25000000
#> 43 State_43  0.23529412
#> 44 State_44  0.30000000
#> 45 State_45  0.22580645
#> 46 State_46  0.33333333
#> 47 State_47  0.29197080
#> 48 State_48  0.33333333
#> 49 State_49  0.32352941
#> 50 State_50  0.18181818
```
