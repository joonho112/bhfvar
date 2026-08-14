.bhf_get_rstan_auto_write <- function() {
  rstan::rstan_options("auto_write")
}

.bhf_set_rstan_auto_write <- function(value) {
  .bhf_assert_flag(value, "auto_write")
  rstan::rstan_options(auto_write = value)
  invisible(value)
}

.bhf_rstan_stan_model <- function(...) rstan::stan_model(...)

.bhf_locate_stan_file <- function() {
  installed <- system.file("stan", "bhf_hybrid.stan", package = "bhfvar")
  candidates <- unique(c(
    installed,
    file.path(getwd(), "inst", "stan", "bhf_hybrid.stan")
  ))
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) return("")
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

.bhf_sha256_file <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path) ||
      !file.exists(path)) {
    stop("Cannot hash a missing Stan model file.", call. = FALSE)
  }
  sha256sum <- get0(
    "sha256sum",
    envir = asNamespace("tools"),
    mode = "function",
    inherits = FALSE
  )
  if (!is.null(sha256sum)) {
    hash <- unname(sha256sum(path))
  } else if (requireNamespace("digest", quietly = TRUE)) {
    hash <- digest::digest(path, algo = "sha256", file = TRUE)
  } else {
    stop(
      "SHA-256 support requires this R version's tools::sha256sum() or ",
      "the 'digest' package.",
      call. = FALSE
    )
  }
  hash <- tolower(as.character(hash))
  if (length(hash) != 1L || !grepl("^[0-9a-f]{64}$", hash)) {
    stop("Failed to compute a valid SHA-256 model fingerprint.", call. = FALSE)
  }
  hash
}

.bhf_model_sha256 <- function() {
  stan_file <- .bhf_locate_stan_file()
  if (!nzchar(stan_file)) {
    stop("Stan model file not found; cannot record model SHA-256.",
         call. = FALSE)
  }
  .bhf_sha256_file(stan_file)
}

#' Compile the BHF Stan Model
#'
#' Compiles the Stan model for the Bayesian Hybrid Framework. This function
#' should be called once at the beginning of each R session before fitting
#' any models.
#'
#' @param verbose Logical. If TRUE, prints compilation progress. Default is TRUE.
#' @param auto_write Logical. If TRUE, saves compiled model to disk for faster
#'   loading in future sessions. Default is TRUE.
#'
#' @return An object of class \code{stanmodel} that can be used with \code{bhf_fit()}.
#'
#' @details
#' This function compiles the Stan model from source code included in the package.
#' Compilation typically takes 30-60 seconds depending on your system.
#'
#' The compiled model is cached by rstan (if \code{auto_write = TRUE}), so
#' subsequent calls to this function in the same or future sessions will be
#' much faster.
#'
#' @section Why Manual Compilation?:
#' Unlike packages that pre-compile Stan models during installation, this package
#' requires explicit compilation. This "defensive" approach:
#' \itemize{
#'   \item Avoids cryptic errors from cached/stale compiled objects
#'   \item Provides clearer error messages when compilation fails
#'   \item Ensures compatibility across different R/Stan/C++ configurations
#'   \item Gives you visibility into what's happening
#' }
#'
#' @section Troubleshooting:
#' If compilation fails:
#' \itemize{
#'   \item Ensure you have a C++ toolchain installed (Rtools on Windows, 
#'     Xcode CLI on macOS, build-essential on Linux)
#'   \item Check that rstan is properly installed: \code{example(stan_model, package = "rstan")}
#'   \item See the package vignette for detailed troubleshooting
#' }
#'
#' @examples
#' \dontrun{
#' # Compile the model (once per session)
#' model <- compile_bhf_model()
#'
#' # Now you can use the model for fitting
#' fit <- bhf_fit(prepared_data, model = model)
#' }
#'
#' @export
compile_bhf_model <- function(verbose = TRUE, auto_write = TRUE) {
  .bhf_assert_flag(verbose, "verbose")
  .bhf_assert_flag(auto_write, "auto_write")

  if (verbose) {
    message("=== Compiling BHF Stan Model ===")
    message("This may take 30-60 seconds on first run...")
  }
  
  old_auto_write <- .bhf_get_rstan_auto_write()
  on.exit(.bhf_set_rstan_auto_write(old_auto_write), add = TRUE)
  .bhf_set_rstan_auto_write(auto_write)

  stan_file <- .bhf_locate_stan_file()
  if (stan_file == "" || !file.exists(stan_file)) {
    stop(
      "Stan model file not found. ",
      "This usually means the package was not installed correctly.\n",
      "Reinstall from an authorized bhfvar source archive or checkout.",
      call. = FALSE
    )
  }
  
  if (verbose) {
    message("Stan file location: ", stan_file)
  }
  
  tryCatch({
    model <- .bhf_rstan_stan_model(
      file = stan_file,
      verbose = verbose
    )
    model_sha256 <- .bhf_sha256_file(stan_file)
    model <- tryCatch({
      attr(model, "bhfvar_model_sha256") <- model_sha256
      attr(model, "bhfvar_model_path") <- stan_file
      model
    }, error = function(e) model)
    
    if (verbose) {
      message("\n=== Compilation Successful ===")
      message("Model is ready for use with bhf_fit()")
    }
    
    return(model)
    
  }, error = function(e) {
    stop(
      "Stan model compilation failed.\n\n",
      "Error message: ", e$message, "\n\n",
      "Common causes:\n",
      "1. Missing C++ toolchain (Rtools on Windows, Xcode CLI on macOS)\n",
      "2. Incompatible rstan/StanHeaders versions\n",
      "3. Insufficient disk space or memory\n\n",
      "Try running: example(stan_model, package = 'rstan')\n",
      "to test if rstan is working correctly.",
      call. = FALSE
    )
  })
}


#' Get Path to Stan Model File
#'
#' Returns the file path to the BHF Stan model included in the package.
#' This is useful for inspecting and hashing the exact bundled model source.
#'
#' @return Character string with the path to the Stan file.
#'
#' @examples
#' \dontrun{
#' # Get the path
#' stan_path <- get_stan_file_path()
#' print(stan_path)
#'
#' # Read and view the model code
#' cat(readLines(stan_path), sep = "\n")
#' }
#'
#' @export
get_stan_file_path <- function() {
  stan_file <- .bhf_locate_stan_file()
  
  if (stan_file == "") {
    stop("Stan model file not found in package.", call. = FALSE)
  }
  
  return(stan_file)
}
