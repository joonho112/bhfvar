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
#'   Higher values (e.g., 0.95) reduce divergent transitions but slow sampling.
#'   Default is 0.90.
#' @param max_treedepth Integer. Maximum tree depth for NUTS. Default is 12.
#' @param refresh Integer. How often to print progress (iterations). Set to 0
#'   to suppress output.
#' @param ... Additional arguments passed to \code{rstan::sampling()}.
#'
#' @return An object of class \code{bhf_fit} containing:
#'   \item{stanfit}{The stanfit object from rstan}
#'   \item{data}{The bhf_data object used for fitting}
#'   \item{model}{The compiled Stan model}
#'   \item{call}{The function call}
#'   \item{diagnostics}{MCMC diagnostics summary}
#'
#' @details
#' The model uses a non-centered parameterization for random effects, which
#' typically improves sampling efficiency. Key parameters are:
#'
#' \describe{
#'   \item{alpha}{Global intercept on logit scale}
#'   \item{sigma_state}{Between-domain standard deviation}
#'   \item{sigma_psu}{Between-PSU standard deviation (within stratum)}
#' }
#'
#' @section Convergence Diagnostics:
#' The function automatically checks for:
#' \itemize{
#'   \item Divergent transitions (should be 0)
#'   \item Low Rhat values (should be < 1.05)
#'   \item Adequate effective sample size (should be > 100)
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
                    adapt_delta = 0.90,
                    max_treedepth = 12,
                    refresh = 200,
                    ...) {
  
  # ==========================================================================
  # Input Validation
  # ==========================================================================
  
  if (!inherits(data, "bhf_data")) {
    stop(
      "'data' must be a 'bhf_data' object from prepare_bhf_data().\n",
      "Got object of class: ", paste(class(data), collapse = ", "),
      call. = FALSE
    )
  }
  
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
  
  # Set cores
  if (is.null(cores)) {
    cores <- min(chains, parallel::detectCores())
  }
  
  # ==========================================================================
  # Sampling
  # ==========================================================================
  
  message("=== Starting BHF Model Fitting ===")
  message("Chains: ", chains, ", Iterations: ", iter, " (warmup: ", warmup, ")")
  message("Observations: ", data$stan_data$N, ", Domains: ", data$stan_data$S)
  
  # Set options
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = cores)
  
  # Run sampling
  stanfit <- rstan::sampling(
    object = model,
    data = data$stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    seed = seed,
    cores = cores,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    refresh = refresh,
    ...
  )
  
  # ==========================================================================
  # Diagnostics
  # ==========================================================================
  
  diagnostics <- compute_diagnostics(stanfit)
  
  # Issue warnings if problems detected
  if (diagnostics$n_divergent > 0) {
    warning(
      diagnostics$n_divergent, " divergent transitions detected.\n",
      "Consider increasing adapt_delta (current: ", adapt_delta, ").",
      call. = FALSE
    )
  }
  
  if (any(diagnostics$rhat_max > 1.05)) {
    warning(
      "Some Rhat values > 1.05, indicating potential convergence issues.\n",
      "Consider increasing iterations.",
      call. = FALSE
    )
  }
  
  if (any(diagnostics$ess_bulk_min < 100)) {
    warning(
      "Low effective sample size detected for some parameters.\n",
      "Consider increasing iterations.",
      call. = FALSE
    )
  }
  
  # ==========================================================================
  # Create Return Object
  # ==========================================================================
  
  result <- list(
    stanfit = stanfit,
    data = data,
    model = model,
    call = match.call(),
    diagnostics = diagnostics
  )
  
  class(result) <- c("bhf_fit", "list")
  
  message("\n=== Model Fitting Complete ===")
  message("Use variance_decomposition() or domain_estimates() to extract results.")
  
  return(result)
}


#' Compute MCMC Diagnostics
#'
#' Internal function to compute diagnostic summaries from a stanfit object.
#'
#' @param stanfit A stanfit object.
#'
#' @return List of diagnostic summaries.
#'
#' @keywords internal
compute_diagnostics <- function(stanfit) {
  
  # Get sampler parameters
  sampler_params <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)
  
  # Divergences
  n_divergent <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
  
  # Tree depth
  max_depth_hits <- sum(sapply(sampler_params, function(x) {
    sum(x[, "treedepth__"] >= 12)  # Assuming default max_treedepth
  }))
  
  # Summary statistics
  fit_summary <- rstan::summary(stanfit)$summary
  
  # Rhat
  rhat_vals <- fit_summary[, "Rhat"]
  rhat_vals <- rhat_vals[!is.na(rhat_vals)]
  
  # Effective sample size
  ess_bulk <- fit_summary[, "n_eff"]
  ess_bulk <- ess_bulk[!is.na(ess_bulk)]
  
  list(
    n_divergent = n_divergent,
    max_depth_hits = max_depth_hits,
    rhat_max = max(rhat_vals, na.rm = TRUE),
    rhat_all_ok = all(rhat_vals < 1.05, na.rm = TRUE),
    ess_bulk_min = min(ess_bulk, na.rm = TRUE),
    ess_bulk_median = median(ess_bulk, na.rm = TRUE)
  )
}


#' @export
print.bhf_fit <- function(x, ...) {
  cat("=== BHF Model Fit ===\n\n")
  
  # Data summary
  cat("Data:\n")
  cat("  Observations: ", x$data$stan_data$N, "\n")
  cat("  Domains: ", x$data$stan_data$S, "\n")
  cat("  Strata: ", x$data$stan_data$H, "\n")
  cat("  PSUs: ", x$data$stan_data$J, "\n")
  
  # Sampling info
  cat("\nSampling:\n")
  cat("  Chains: ", x$stanfit@sim$chains, "\n")
  cat("  Iterations: ", x$stanfit@sim$iter, "\n")
  cat("  Warmup: ", x$stanfit@sim$warmup, "\n")
  
  # Diagnostics
  cat("\nDiagnostics:\n")
  cat("  Divergent transitions: ", x$diagnostics$n_divergent, "\n")
  cat("  Max Rhat: ", sprintf("%.3f", x$diagnostics$rhat_max), "\n")
  cat("  Min ESS: ", sprintf("%.0f", x$diagnostics$ess_bulk_min), "\n")
  
  # Key parameters
  cat("\nKey Parameters (posterior mean [95% CI]):\n")
  summ <- rstan::summary(x$stanfit, pars = c("alpha", "sigma_state", "sigma_psu"))$summary
  for (param in c("alpha", "sigma_state", "sigma_psu")) {
    if (param %in% rownames(summ)) {
      cat(sprintf("  %s: %.3f [%.3f, %.3f]\n",
                  param,
                  summ[param, "mean"],
                  summ[param, "2.5%"],
                  summ[param, "97.5%"]))
    }
  }
  
  cat("\nUse variance_decomposition(fit) or domain_estimates(fit) for full results.\n")
  
  invisible(x)
}


#' @export
summary.bhf_fit <- function(object, ...) {
  
  cat("=== BHF Model Summary ===\n\n")
  
  # Extract key quantities
  vd <- variance_decomposition(object, print = FALSE)
  
  cat("--- Variance Decomposition ---\n\n")
  
  cat("Logit Scale (Estimand A - Policy):\n")
  cat(sprintf("  Between-domain variance: %.4f\n", vd$logit$var_between_mean))
  cat(sprintf("  Within-domain variance: %.4f\n", vd$logit$var_within_mean))
  cat(sprintf("  ICC: %.3f [%.3f, %.3f]\n", 
              vd$logit$icc_mean,
              vd$logit$icc_q025,
              vd$logit$icc_q975))
  
  cat("\nProbability Scale (Estimand B - Descriptive):\n")
  cat(sprintf("  Between-domain variance: %.4f\n", vd$prob$var_between_mean))
  cat(sprintf("  Within-domain variance: %.4f\n", vd$prob$var_within_mean))
  cat(sprintf("  ICC: %.3f [%.3f, %.3f]\n", 
              vd$prob$icc_mean,
              vd$prob$icc_q025,
              vd$prob$icc_q975))
  
  if (object$data$stan_data$use_deattenuation == 1) {
    cat("\nDe-attenuated (Estimand A* - Policy Adjusted):\n")
    cat(sprintf("  Between-domain variance: %.4f\n", vd$deatten$var_between_mean))
    cat(sprintf("  ICC: %.3f [%.3f, %.3f]\n",
                vd$deatten$icc_mean,
                vd$deatten$icc_q025,
                vd$deatten$icc_q975))
  }
  
  cat("\n--- Model Diagnostics ---\n\n")
  cat(sprintf("  Divergent transitions: %d\n", object$diagnostics$n_divergent))
  cat(sprintf("  Max Rhat: %.3f\n", object$diagnostics$rhat_max))
  cat(sprintf("  Min ESS (bulk): %.0f\n", object$diagnostics$ess_bulk_min))
  
  invisible(object)
}
