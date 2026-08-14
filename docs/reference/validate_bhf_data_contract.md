# Validate the Complete Prepared-Data Contract

Internal final gate for cross-field transport, mapping, and provenance
invariants. Field-local Stan support is delegated to
[`validate_stan_data()`](https://joonho112.github.io/bhfvar/reference/validate_stan_data.md).

## Usage

``` r
validate_bhf_data_contract(x)
```

## Arguments

- x:

  A prepared `bhf_data` object.

## Value

`x`, invisibly, or a classed `bhf_data_contract_error`.
