# Get Path to Stan Model File

Returns the file path to the BHF Stan model included in the package.
This is useful if you want to inspect or modify the model code.

## Usage

``` r
get_stan_file_path()
```

## Value

Character string with the path to the Stan file.

## Examples

``` r
# Get the path
stan_path <- get_stan_file_path()
print(stan_path)
#> [1] "/private/var/folders/22/41f2lr_j0rz4yfdj76vg19m80000gq/T/Rtmp5jBq8z/temp_libpath451b2e9be6fa/bhfvar/stan/bhf_hybrid.stan"

# Read and view the model code
if (FALSE) { # \dontrun{
cat(readLines(stan_path), sep = "\n")
} # }
```
