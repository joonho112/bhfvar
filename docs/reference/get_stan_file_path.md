# Get Path to Stan Model File

Returns the file path to the BHF Stan model included in the package.
This is useful for inspecting and hashing the exact bundled model
source.

## Usage

``` r
get_stan_file_path()
```

## Value

Character string with the path to the Stan file.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get the path
stan_path <- get_stan_file_path()
print(stan_path)

# Read and view the model code
cat(readLines(stan_path), sep = "\n")
} # }
```
