# Compile the BHF Stan Model

Compiles the Stan model for the Bayesian Hybrid Framework. This function
should be called once at the beginning of each R session before fitting
any models.

## Usage

``` r
compile_bhf_model(verbose = TRUE, auto_write = TRUE)
```

## Arguments

- verbose:

  Logical. If TRUE, prints compilation progress. Default is TRUE.

- auto_write:

  Logical. If TRUE, saves compiled model to disk for faster loading in
  future sessions. Default is TRUE.

## Value

An object of class `stanmodel` that can be used with
[`bhf_fit()`](https://joonho112.github.io/bhfvar/reference/bhf_fit.md).

## Details

This function compiles the Stan model from source code included in the
package. Compilation typically takes 30-60 seconds depending on your
system.

The compiled model is cached by rstan (if `auto_write = TRUE`), so
subsequent calls to this function in the same or future sessions will be
much faster.

## Why Manual Compilation?

Unlike packages that pre-compile Stan models during installation, this
package requires explicit compilation. This "defensive" approach:

- Avoids cryptic errors from cached/stale compiled objects

- Provides clearer error messages when compilation fails

- Ensures compatibility across different R/Stan/C++ configurations

- Gives you visibility into what's happening

## Troubleshooting

If compilation fails:

- Ensure you have a C++ toolchain installed (Rtools on Windows, Xcode
  CLI on macOS, build-essential on Linux)

- Check that rstan is properly installed:
  `example(stan_model, package = "rstan")`

- See the package vignette for detailed troubleshooting

## Examples

``` r
if (FALSE) { # \dontrun{
# Compile the model (once per session)
model <- compile_bhf_model()

# Now you can use the model for fitting
fit <- bhf_fit(prepared_data, model = model)
} # }
```
