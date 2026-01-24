# Validate Input Data

Internal function to validate input data for prepare_bhf_data().

## Usage

``` r
validate_input_data(data, outcome, domain, strata, psu, weights)
```

## Arguments

- data:

  Data frame to validate.

- outcome:

  Name of outcome variable.

- domain:

  Name of domain variable.

- strata:

  Name of strata variable.

- psu:

  Name of PSU variable.

- weights:

  Name of weight variable.

## Value

NULL (invisibly). Throws an error if validation fails.
