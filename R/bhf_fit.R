.bhf_fit_contract_spec <- function() {
  list(
    schema_version = "0.5.0",
    contract_id = "bhfvar-fit-contract-0.5.0",
    package_version = .bhf_namespace_version("bhfvar")
  )
}

.bhf_assert_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be one non-missing logical value.", call. = FALSE)
  }
  invisible(x)
}

.bhf_assert_integer_scalar <- function(x, name, lower = 0L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != as.integer(x) || x < lower || x > .Machine$integer.max) {
    stop("`", name, "` must be one integer >= ", lower, ".", call. = FALSE)
  }
  as.integer(x)
}

.bhf_assert_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x <= 0 || x >= 1) {
    stop("`", name, "` must be one finite number strictly between 0 and 1.",
         call. = FALSE)
  }
  as.numeric(x)
}

.bhf_validate_sampling_dots <- function(dots) {
  if (!length(dots)) return(character())
  dot_names <- names(dots)
  if (is.null(dot_names) || anyNA(dot_names) || any(!nzchar(dot_names))) {
    stop("All additional sampling arguments in `...` must be named.",
         call. = FALSE)
  }
  if (anyDuplicated(dot_names)) {
    stop("Additional sampling arguments in `...` must have unique names.",
         call. = FALSE)
  }
  reserved <- c(
    "object", "data", "chains", "iter", "warmup", "seed", "cores",
    "control", "refresh"
  )
  duplicate <- intersect(dot_names, reserved)
  if (length(duplicate)) {
    stop(
      "Do not pass reserved sampling argument(s) through `...`: ",
      paste(duplicate, collapse = ", "),
      ". Use the corresponding `bhf_fit()` formal argument; the Stan ",
      "`control` list is built from `adapt_delta` and `max_treedepth`.",
      call. = FALSE
    )
  }
  dot_names
}

.bhf_rstan_sampling <- function(...) rstan::sampling(...)
.bhf_rstan_summary <- function(...) rstan::summary(...)
.bhf_rstan_sampler_params <- function(...) rstan::get_sampler_params(...)

.bhf_mc_cores_state <- function() {
  list(
    present = "mc.cores" %in% names(options()),
    value = getOption("mc.cores")
  )
}

.bhf_restore_mc_cores <- function(state) {
  if (isTRUE(state$present)) {
    options(mc.cores = state$value)
  } else {
    options(mc.cores = NULL)
  }
  invisible(NULL)
}

.bhf_namespace_version <- function(package) {
  tryCatch(
    as.character(getNamespaceVersion(package)),
    error = function(e) NA_character_
  )
}

.bhf_session_metadata <- function(started_at, completed_at = Sys.time()) {
  sys <- Sys.info()
  list(
    started_at_utc = format(started_at, tz = "UTC", usetz = TRUE),
    completed_at_utc = format(completed_at, tz = "UTC", usetz = TRUE),
    r_version = R.version.string,
    platform = R.version$platform,
    os = unname(sys[["sysname"]]),
    release = unname(sys[["release"]]),
    bhfvar_version = .bhf_namespace_version("bhfvar"),
    rstan_version = .bhf_namespace_version("rstan")
  )
}

.bhf_fit_provenance <- function(data, model_sha256, sampling, started_at) {
  dimensions <- data$stan_data[c("N", "S", "H", "J")]
  list(
    contract = c(
      .bhf_fit_contract_spec(),
      list(
        data_schema_version = data$schema_version,
        data_contract_id = data$contract_id,
        stan_data_schema_version = data$stan_data_schema_version,
        model_sha256 = model_sha256
      )
    ),
    data = list(
      schema_version = data$schema_version,
      contract_id = data$contract_id,
      stan_data_schema_version = data$stan_data_schema_version,
      dimensions = dimensions,
      n_original = data$row_provenance$n_original,
      n_used = data$row_provenance$n_used,
      retained_rows = data$row_provenance$retained_rows,
      dropped_rows = data$row_provenance$dropped_rows
    ),
    weights = c(
      data$provenance$weights,
      list(
        original_row = data$weight_info$original_row,
        raw = data$weight_info$raw,
        w_lik = data$weight_info$likelihood
      )
    ),
    population_shares = c(
      data$provenance$population_shares,
      list(values = data$population_share_info$values)
    ),
    sampling_variances = c(
      data$provenance$sampling_variances,
      list(
        enabled = data$sampling_variance_info$enabled,
        named_values = data$sampling_variance_info$named_values,
        stan_values = data$sampling_variance_info$stan_values
      )
    ),
    prior = data$provenance$prior,
    sampling = sampling,
    session = .bhf_session_metadata(started_at)
  )
}

#' Fit the Bayesian Hybrid Framework Model
#'
#' Fits the BHF model using Hamiltonian Monte Carlo via Stan.
#'
#' @param data An object of class \code{bhf_data} from \code{prepare_bhf_data()}.
#' @param model A compiled Stan model from \code{compile_bhf_model()}. If NULL,
#'   the function will attempt to compile the model (with a warning).
#' @param chains Integer. Number of MCMC chains. Default is 4.
#' @param iter Integer. Total number of iterations per chain (including warmup).
#'   Default is 2000.
#' @param warmup Integer. Number of warmup iterations per chain. Default is
#'   half of \code{iter}.
#' @param seed Integer. Random seed for reproducibility.
#' @param cores Integer. Number of cores for parallel chains. Default is
#'   \code{min(chains, parallel::detectCores())}.
#' @param adapt_delta Numeric between 0 and 1. Target acceptance probability.
#'   Higher values reduce divergent transitions but slow sampling. Default is
#'   0.95.
#' @param max_treedepth Integer. Maximum tree depth for NUTS. Default is 12.
#' @param refresh Integer. How often to print progress (iterations). Set to 0
#'   to suppress output.
#' @param ... Additional arguments passed to \code{rstan::sampling()}.
#'
#' @return An object of class \code{bhf_fit} containing:
#'   \item{schema_version}{Version of the fit-object contract}
#'   \item{contract_id}{Stable fit-object contract identifier}
#'   \item{package_version}{Package version that created the fit}
#'   \item{model_sha256}{SHA-256 fingerprint of the Stan source}
#'   \item{stanfit}{The stanfit object from rstan}
#'   \item{data}{The bhf_data object used for fitting}
#'   \item{model}{The compiled Stan model}
#'   \item{call}{The function call}
#'   \item{diagnostics}{MCMC diagnostics summary}
#'   \item{provenance}{Data, weight, population-share, sampling-variance,
#'     prior, sampling-control, and session metadata}
#'
#' @details
#' The model uses a non-centered parameterization for random effects, which
#' typically improves sampling efficiency. Key parameters are:
#'
#' \describe{
#'   \item{alpha}{Global intercept on logit scale}
#'   \item{sigma_state}{Between-domain standard deviation}
#'   \item{sigma_stratum}{Between-stratum design standard deviation}
#'   \item{sigma_psu}{Between-PSU standard deviation (within stratum)}
#' }
#'
#' @section Convergence Diagnostics:
#' The function automatically checks for:
#' \itemize{
#'   \item Divergent transitions (should be 0)
#'   \item Low Rhat values (should be < 1.05)
#'   \item Legacy rstan \code{n_eff} (should be > 100; this is not labeled as
#'     bulk or tail ESS)
#'   \item Tree depth saturation
#' }
#' Warnings are issued if any diagnostics suggest convergence problems.
#'
#' @section Typical Runtime:
#' With default settings (4 chains, 2000 iterations):
#' \itemize{
#'   \item Small data (N < 500): 1-2 minutes
#'   \item Medium data (N ~ 3000): 5-10 minutes
#'   \item Large data (N > 10000): 20+ minutes
#' }
#'
#' @examples
#' \dontrun{
#' # Step 1: Compile model (once per session)
#' model <- compile_bhf_model()
#'
#' # Step 2: Prepare data
#' data(bhf_synthetic_data)
#' prepared <- prepare_bhf_data(
#'   bhf_synthetic_data,
#'   outcome = "has_subsidy",
#'   domain = "state",
#'   strata = "stratum",
#'   psu = "psu",
#'   weights = "weight"
#' )
#'
#' # Step 3: Fit model
#' fit <- bhf_fit(prepared, model = model, chains = 4, iter = 2000)
#'
#' # Step 4: Examine results
#' print(fit)
#' summary(fit)
#' }
#'
#' @export
bhf_fit <- function(data,
                    model = NULL,
                    chains = 4,
                    iter = 2000,
                    warmup = floor(iter / 2),
                    seed = 1234,
                    cores = NULL,
                    adapt_delta = 0.95,
                    max_treedepth = 12,
                    refresh = 200,
                    ...) {
  started_at <- Sys.time()
  dots <- list(...)
  dot_names <- .bhf_validate_sampling_dots(dots)

  if (!inherits(data, "bhf_data")) {
    stop(
      "'data' must be a 'bhf_data' object from prepare_bhf_data().\n",
      "Got object of class: ", paste(class(data), collapse = ", "),
      call. = FALSE
    )
  }
  if (exists(".bhf_legacy_assert_current_bhf_data", mode = "function")) {
    .bhf_legacy_assert_current_bhf_data(data)
  } else {
    validate_bhf_data_contract(data)
  }

  chains <- .bhf_assert_integer_scalar(chains, "chains", 1L)
  iter <- .bhf_assert_integer_scalar(iter, "iter", 1L)
  warmup <- .bhf_assert_integer_scalar(warmup, "warmup", 0L)
  seed <- .bhf_assert_integer_scalar(seed, "seed", 1L)
  max_treedepth <- .bhf_assert_integer_scalar(
    max_treedepth, "max_treedepth", 1L
  )
  refresh <- .bhf_assert_integer_scalar(refresh, "refresh", 0L)
  adapt_delta <- .bhf_assert_probability(adapt_delta, "adapt_delta")
  if (warmup >= iter) {
    stop("`warmup` must be strictly less than `iter`.", call. = FALSE)
  }
  if (is.null(cores)) {
    detected <- suppressWarnings(parallel::detectCores(logical = FALSE))
    if (length(detected) != 1L || is.na(detected) || detected < 1L) detected <- 1L
    cores <- min(chains, as.integer(detected))
  }
  cores <- .bhf_assert_integer_scalar(cores, "cores", 1L)

  if (is.null(model)) {
    warning(
      "No compiled model provided. Compiling model now...\n",
      "For faster execution, compile once with compile_bhf_model() and reuse.",
      call. = FALSE
    )
    model <- compile_bhf_model(verbose = FALSE)
  }
  
  if (!inherits(model, "stanmodel")) {
    stop(
      "'model' must be a compiled Stan model from compile_bhf_model().\n",
      "Got object of class: ", paste(class(model), collapse = ", "),
      call. = FALSE
    )
  }

  model_sha256 <- attr(model, "bhfvar_model_sha256", exact = TRUE)
  if (!is.character(model_sha256) || length(model_sha256) != 1L ||
      !grepl("^[0-9a-f]{64}$", model_sha256)) {
    model_sha256 <- .bhf_model_sha256()
  }
  control <- list(
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )

  message("=== Starting BHF Model Fitting ===")
  message("Chains: ", chains, ", Iterations: ", iter, " (warmup: ", warmup, ")")
  message("Observations: ", data$stan_data$N, ", Domains: ", data$stan_data$S)

  old_auto_write <- .bhf_get_rstan_auto_write()
  old_mc_cores <- .bhf_mc_cores_state()
  on.exit({
    .bhf_set_rstan_auto_write(old_auto_write)
    .bhf_restore_mc_cores(old_mc_cores)
  }, add = TRUE)
  .bhf_set_rstan_auto_write(TRUE)
  options(mc.cores = cores)

  sampling_args <- c(
    list(
      object = model,
      data = data$stan_data,
      chains = chains,
      iter = iter,
      warmup = warmup,
      seed = seed,
      cores = cores,
      control = control,
      refresh = refresh
    ),
    dots
  )
  stanfit <- do.call(.bhf_rstan_sampling, sampling_args)

  structural_parameters <- c(
    "alpha", "sigma_state", "sigma_stratum", "sigma_psu"
  )
  diagnostics <- compute_diagnostics(
    stanfit,
    max_treedepth = max_treedepth,
    pars = structural_parameters
  )

  if (diagnostics$n_divergent > 0) {
    warning(
      diagnostics$n_divergent, " divergent transitions detected.\n",
      "Consider increasing adapt_delta (current: ", adapt_delta, ").",
      call. = FALSE
    )
  }
  if (isTRUE(diagnostics$rhat_max >= 1.05)) {
    warning(
      "Some Rhat values >= 1.05, indicating potential convergence issues.\n",
      "Consider increasing iterations.",
      call. = FALSE
    )
  }
  if (isTRUE(diagnostics$n_eff_min < 100)) {
    warning(
      "Low legacy rstan n_eff detected for some structural parameters.\n",
      "Consider increasing iterations.",
      call. = FALSE
    )
  }
  if (diagnostics$max_depth_hits > 0) {
    warning(
      diagnostics$max_depth_hits,
      " post-warmup transitions reached the configured max_treedepth (",
      max_treedepth, ").",
      call. = FALSE
    )
  }

  sampling_provenance <- list(
    seed = seed,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    refresh = refresh,
    control = control,
    additional_argument_names = dot_names,
    additional_arguments = dots
  )
  provenance <- .bhf_fit_provenance(
    data = data,
    model_sha256 = model_sha256,
    sampling = sampling_provenance,
    started_at = started_at
  )
  contract <- .bhf_fit_contract_spec()
  result <- list(
    schema_version = contract$schema_version,
    contract_id = contract$contract_id,
    package_version = contract$package_version,
    model_sha256 = model_sha256,
    stanfit = stanfit,
    data = data,
    model = model,
    call = match.call(),
    diagnostics = diagnostics,
    provenance = provenance
  )
  class(result) <- c("bhf_fit", "list")

  message("\n=== Model Fitting Complete ===")
  message("Use variance_decomposition() or domain_estimates() to extract results.")
  
  result
}


#' Compute MCMC Diagnostics
#'
#' Internal function to compute diagnostic summaries from a stanfit object.
#'
#' @param stanfit A stanfit object.
#' @param max_treedepth Configured maximum tree depth used for this fit.
#' @param pars Structural parameter names to diagnose. Generated quantities are
#'   intentionally excluded.
#'
#' @return List of diagnostic summaries.
#'
#' @keywords internal
compute_diagnostics <- function(
    stanfit,
    max_treedepth = 12L,
    pars = c("alpha", "sigma_state", "sigma_stratum", "sigma_psu")) {
  max_treedepth <- .bhf_assert_integer_scalar(
    max_treedepth, "max_treedepth", 1L
  )
  if (!is.character(pars) || !length(pars) || anyNA(pars) ||
      any(!nzchar(pars)) || anyDuplicated(pars)) {
    stop("`pars` must be a non-empty vector of unique parameter names.",
         call. = FALSE)
  }

  sampler_params <- .bhf_rstan_sampler_params(
    stanfit,
    inc_warmup = FALSE
  )
  if (!is.list(sampler_params)) sampler_params <- list()
  sampler_column <- function(name) {
    values <- unlist(lapply(sampler_params, function(chain) {
      if (!is.matrix(chain) || !name %in% colnames(chain)) return(numeric())
      as.numeric(chain[, name])
    }), use.names = FALSE)
    values[is.finite(values)]
  }
  divergent <- sampler_column("divergent__")
  treedepth <- sampler_column("treedepth__")

  summary_result <- .bhf_rstan_summary(stanfit, pars = pars)
  fit_summary <- summary_result$summary
  if (!is.matrix(fit_summary)) {
    fit_summary <- matrix(numeric(), nrow = 0L, ncol = 0L)
  }
  if (!is.null(rownames(fit_summary))) {
    keep <- intersect(pars, rownames(fit_summary))
    fit_summary <- fit_summary[keep, , drop = FALSE]
  } else {
    fit_summary <- fit_summary[0, , drop = FALSE]
  }
  finite_summary_column <- function(name) {
    if (!name %in% colnames(fit_summary)) return(numeric())
    values <- as.numeric(fit_summary[, name])
    values[is.finite(values)]
  }
  rhat <- finite_summary_column("Rhat")
  n_eff <- finite_summary_column("n_eff")
  finite_or_missing <- function(values, fun) {
    if (!length(values)) NA_real_ else as.numeric(fun(values))
  }

  list(
    parameters = rownames(fit_summary),
    n_divergent = as.integer(sum(divergent > 0)),
    max_treedepth_configured = max_treedepth,
    max_treedepth_observed = if (length(treedepth)) {
      as.integer(max(treedepth))
    } else {
      NA_integer_
    },
    max_depth_hits = as.integer(sum(treedepth >= max_treedepth)),
    rhat_max = finite_or_missing(rhat, max),
    rhat_all_ok = if (length(rhat)) all(rhat < 1.05) else NA,
    rhat_n = length(rhat),
    n_eff_min = finite_or_missing(n_eff, min),
    n_eff_median = finite_or_missing(n_eff, stats::median),
    n_eff_n = length(n_eff),
    n_eff_label = "legacy rstan n_eff (not bulk or tail ESS)"
  )
}


#' @export
print.bhf_fit <- function(x, ...) {
  .bhf_validate_result_fit(x)
  sampling <- x$provenance$sampling
  diagnostics <- x$diagnostics
  value_or_na <- function(value, digits = 3) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      return("NA")
    }
    formatC(value, format = "f", digits = digits)
  }

  cat("=== BHF Model Fit (schema ", x$schema_version, ") ===\n\n",
      sep = "")
  cat("Data: ", x$data$stan_data$N, " observations; ",
      x$data$stan_data$S, " domains; ", x$data$stan_data$H,
      " strata; ", x$data$stan_data$J, " PSUs\n", sep = "")
  cat("Model SHA-256: ", x$model_sha256, "\n", sep = "")
  cat("Sampling: ", sampling$chains, " chains x ", sampling$iter,
      " iterations (", sampling$warmup, " warmup); seed ", sampling$seed,
      "\n", sep = "")
  cat("Control: adapt_delta=", sampling$control$adapt_delta,
      ", max_treedepth=", sampling$control$max_treedepth, "\n", sep = "")
  cat("Diagnostics: divergences=", diagnostics$n_divergent,
      ", treedepth hits=", diagnostics$max_depth_hits,
      ", max Rhat=", value_or_na(diagnostics$rhat_max),
      ", min legacy n_eff=", value_or_na(diagnostics$n_eff_min, 0),
      "\n", sep = "")
  cat("A*: ", if (x$data$stan_data$use_deattenuation == 1L) {
    "available from fixed sampling variances"
  } else "unavailable (deattenuation='none')", "\n", sep = "")
  cat("Use variance_decomposition(), domain_estimates(), or overall_estimate().\n")
  invisible(x)
}


#' @export
summary.bhf_fit <- function(object, ...) {
  .bhf_validate_result_fit(object)
  cat("=== BHF Model Summary ===\n\n")
  print.bhf_fit(object)
  cat("\n--- Article-aligned variance results ---\n")
  decomposition <- variance_decomposition(object, print = FALSE)
  print(decomposition)
  invisible(object)
}
