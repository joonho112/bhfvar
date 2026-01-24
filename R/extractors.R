#' Extract Variance Decomposition Results
#'
#' Extracts and summarizes the variance decomposition from a fitted BHF model.
#' Provides estimates for multiple estimands: Policy (A), Descriptive (B),
#' and De-attenuated (A*).
#'
#' @param fit An object of class \code{bhf_fit} from \code{bhf_fit()}.
#' @param prob Numeric. Probability for credible intervals. Default is 0.95.
#' @param print Logical. If TRUE, prints a formatted summary. Default is TRUE.
#'
#' @return A list with components:
#'   \item{logit}{Variance components on logit scale (Estimand A)}
#'   \item{prob}{Variance components on probability scale (Estimand B)}
#'   \item{deatten}{De-attenuated variance components (Estimand A*)}
#'   \item{summary_table}{Formatted summary table}
#'
#' @details
#' The function extracts three sets of variance components:
#'
#' \describe{
#'   \item{Estimand A (Policy, logit scale)}{
#'     Variance decomposition on the latent logit scale.
#'     \code{icc_state} = sigma_state^2 / (sigma_state^2 + sigma_psu^2 + pi^2/3)
#'   }
#'   \item{Estimand B (Descriptive, probability scale)}{
#'     Variance decomposition on the observed probability scale.
#'     \code{icc_prob = Var(p_s) / (Var(p_s) + E(p_s(1-p_s)))}
#'   }
#'   \item{Estimand A* (Policy adjusted, de-attenuated)}{
#'     Estimand B with finite-sample variance inflation removed.
#'     \code{icc_deatten = (Var(p_s) - V_hat) / (Var(p_s) - V_hat + E(p_s(1-p_s)))}
#'   }
#' }
#'
#' @section Interpretation:
#' The ICC (Intraclass Correlation Coefficient) represents the proportion
#' of total variance attributable to between-domain differences. Higher
#' values indicate more geographic heterogeneity.
#'
#' The difference between ICC_B and ICC_A* represents the amount of
#' apparent heterogeneity that is actually due to sampling noise rather
#' than true substantive differences.
#'
#' @examples
#' \dontrun{
#' # After fitting
#' fit <- bhf_fit(prepared_data, model = model)
#'
#' # Get variance decomposition
#' vd <- variance_decomposition(fit)
#'
#' # Access specific components
#' vd$logit$icc_mean      # ICC on logit scale
#' vd$prob$icc_mean       # ICC on probability scale
#' vd$deatten$icc_mean    # De-attenuated ICC
#' }
#'
#' @export
variance_decomposition <- function(fit, prob = 0.95, print = TRUE) {
  
  if (!inherits(fit, "bhf_fit")) {
    stop("'fit' must be a 'bhf_fit' object from bhf_fit()", call. = FALSE)
  }
  
  # Quantile probabilities
  alpha <- (1 - prob) / 2
  probs <- c(alpha, 0.5, 1 - alpha)
  
  # Extract posterior samples
  samples <- rstan::extract(fit$stanfit)
  
  # ========================================================================
  # Logit Scale (Estimand A)
  # ========================================================================
  
  var_between_logit <- samples$var_between_state
  var_within_logit <- samples$var_within_state
  icc_logit <- samples$icc_state
  
  logit_results <- list(
    var_between_mean = mean(var_between_logit),
    var_between_q025 = quantile(var_between_logit, probs[1]),
    var_between_q975 = quantile(var_between_logit, probs[3]),
    var_within_mean = mean(var_within_logit),
    var_within_q025 = quantile(var_within_logit, probs[1]),
    var_within_q975 = quantile(var_within_logit, probs[3]),
    icc_mean = mean(icc_logit),
    icc_sd = sd(icc_logit),
    icc_q025 = quantile(icc_logit, probs[1]),
    icc_q500 = quantile(icc_logit, probs[2]),
    icc_q975 = quantile(icc_logit, probs[3])
  )
  
  # ========================================================================
  # Probability Scale (Estimand B)
  # ========================================================================
  
  var_between_prob <- samples$var_between_prob
  var_within_prob <- samples$var_within_prob
  icc_prob <- samples$icc_prob
  
  prob_results <- list(
    var_between_mean = mean(var_between_prob),
    var_between_q025 = quantile(var_between_prob, probs[1]),
    var_between_q975 = quantile(var_between_prob, probs[3]),
    var_within_mean = mean(var_within_prob),
    var_within_q025 = quantile(var_within_prob, probs[1]),
    var_within_q975 = quantile(var_within_prob, probs[3]),
    icc_mean = mean(icc_prob),
    icc_sd = sd(icc_prob),
    icc_q025 = quantile(icc_prob, probs[1]),
    icc_q500 = quantile(icc_prob, probs[2]),
    icc_q975 = quantile(icc_prob, probs[3])
  )
  
  # ========================================================================
  # De-attenuated (Estimand A*)
  # ========================================================================
  
  var_between_deatten <- samples$var_between_deatten
  icc_deatten <- samples$icc_deatten
  
  deatten_results <- list(
    var_between_mean = mean(var_between_deatten),
    var_between_q025 = quantile(var_between_deatten, probs[1]),
    var_between_q975 = quantile(var_between_deatten, probs[3]),
    icc_mean = mean(icc_deatten),
    icc_sd = sd(icc_deatten),
    icc_q025 = quantile(icc_deatten, probs[1]),
    icc_q500 = quantile(icc_deatten, probs[2]),
    icc_q975 = quantile(icc_deatten, probs[3])
  )
  
  # ========================================================================
  # Summary Table
  # ========================================================================
  
  summary_table <- data.frame(
    Estimand = c("A (Logit)", "B (Probability)", "A* (De-attenuated)"),
    ICC_Mean = c(logit_results$icc_mean, 
                 prob_results$icc_mean, 
                 deatten_results$icc_mean),
    ICC_SD = c(logit_results$icc_sd, 
               prob_results$icc_sd, 
               deatten_results$icc_sd),
    ICC_Lower = c(logit_results$icc_q025, 
                  prob_results$icc_q025, 
                  deatten_results$icc_q025),
    ICC_Upper = c(logit_results$icc_q975, 
                  prob_results$icc_q975, 
                  deatten_results$icc_q975),
    stringsAsFactors = FALSE
  )
  
  # ========================================================================
  # Print if requested
  # ========================================================================
  
  if (print) {
    cat("=== Variance Decomposition Results ===\n\n")
    
    cat("Estimand A (Policy - Logit Scale):\n")
    cat(sprintf("  ICC: %.3f (95%% CI: [%.3f, %.3f])\n",
                logit_results$icc_mean,
                logit_results$icc_q025,
                logit_results$icc_q975))
    cat(sprintf("  Var(between): %.4f, Var(within): %.4f\n\n",
                logit_results$var_between_mean,
                logit_results$var_within_mean))
    
    cat("Estimand B (Descriptive - Probability Scale):\n")
    cat(sprintf("  ICC: %.3f (95%% CI: [%.3f, %.3f])\n",
                prob_results$icc_mean,
                prob_results$icc_q025,
                prob_results$icc_q975))
    cat(sprintf("  Var(between): %.4f, Var(within): %.4f\n\n",
                prob_results$var_between_mean,
                prob_results$var_within_mean))
    
    if (fit$data$stan_data$use_deattenuation == 1) {
      cat("Estimand A* (Policy Adjusted - De-attenuated):\n")
      cat(sprintf("  ICC: %.3f (95%% CI: [%.3f, %.3f])\n",
                  deatten_results$icc_mean,
                  deatten_results$icc_q025,
                  deatten_results$icc_q975))
      cat(sprintf("  Var(between): %.4f\n\n",
                  deatten_results$var_between_mean))
      
      # Interpretation
      reduction <- (prob_results$icc_mean - deatten_results$icc_mean) / 
                   prob_results$icc_mean * 100
      cat("Interpretation:\n")
      cat(sprintf("  %.1f%% of observed between-domain variance is sampling noise\n",
                  reduction))
    }
  }
  
  # ========================================================================
  # Return
  # ========================================================================
  
  result <- list(
    logit = logit_results,
    prob = prob_results,
    deatten = deatten_results,
    summary_table = summary_table
  )
  
  class(result) <- c("bhf_variance_decomposition", "list")
  
  invisible(result)
}


#' Extract Domain-Level Estimates
#'
#' Extracts posterior summaries of domain-specific probabilities from a
#' fitted BHF model.
#'
#' @param fit An object of class \code{bhf_fit} from \code{bhf_fit()}.
#' @param type Character. Type of probability to extract:
#'   \itemize{
#'     \item \code{"marginal"}: Marginal probabilities (integrating out within-domain variation)
#'     \item \code{"conditional"}: Conditional probabilities (given domain random effect)
#'   }
#'   Default is "marginal".
#' @param prob Numeric. Probability for credible intervals. Default is 0.95.
#'
#' @return A data frame with columns:
#'   \item{domain}{Domain label (from original data)}
#'   \item{domain_id}{Domain ID (1:S)}
#'   \item{mean}{Posterior mean}
#'   \item{sd}{Posterior standard deviation}
#'   \item{q025}{Lower credible interval bound}
#'   \item{q500}{Posterior median}
#'   \item{q975}{Upper credible interval bound}
#'   \item{pop_share}{Population share of domain}
#'   \item{reliability}{Reliability/shrinkage factor for domain}
#'
#' @details
#' The two types of probabilities differ in their interpretation:
#'
#' \describe{
#'   \item{Marginal probabilities}{Average probability for a randomly selected
#'     individual from the domain, integrating over all uncertainty. These are
#'     appropriate for population-level inference.}
#'   \item{Conditional probabilities}{Probability given the estimated domain
#'     random effect. These represent the model's "best guess" for the domain
#'     but ignore uncertainty in the random effect.}
#' }
#'
#' @examples
#' \dontrun{
#' # After fitting
#' fit <- bhf_fit(prepared_data, model = model)
#'
#' # Get domain estimates
#' estimates <- domain_estimates(fit, type = "marginal")
#'
#' # View results
#' head(estimates)
#'
#' # Plot estimates
#' library(ggplot2)
#' ggplot(estimates, aes(x = reorder(domain, mean), y = mean)) +
#'   geom_point() +
#'   geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.2) +
#'   coord_flip() +
#'   labs(x = "Domain", y = "Probability")
#' }
#'
#' @export
domain_estimates <- function(fit, type = c("marginal", "conditional"), prob = 0.95) {
  
  if (!inherits(fit, "bhf_fit")) {
    stop("'fit' must be a 'bhf_fit' object from bhf_fit()", call. = FALSE)
  }
  
  type <- match.arg(type)
  
  # Quantile probabilities
  alpha <- (1 - prob) / 2
  probs <- c(alpha, 0.5, 1 - alpha)
  
  # Extract posterior samples
  samples <- rstan::extract(fit$stanfit)
  
  # Select appropriate parameter
  if (type == "marginal") {
    p_samples <- samples$p_state_marginal  # n_iter x S matrix
  } else {
    p_samples <- samples$p_state_conditional
  }
  
  # Get number of domains
  S <- ncol(p_samples)
  
  # Compute summaries for each domain
  results <- data.frame(
    domain_id = 1:S,
    mean = apply(p_samples, 2, mean),
    sd = apply(p_samples, 2, sd),
    q025 = apply(p_samples, 2, quantile, probs = probs[1]),
    q500 = apply(p_samples, 2, quantile, probs = probs[2]),
    q975 = apply(p_samples, 2, quantile, probs = probs[3]),
    stringsAsFactors = FALSE
  )
  
  # Add domain labels
  domain_mapping <- fit$data$mapping$domain
  results$domain <- domain_mapping$label[match(results$domain_id, domain_mapping$id)]
  
  # Add population shares
  results$pop_share <- fit$data$stan_data$w_state_pop_share
  
  # Add reliability
  if ("reliability_state" %in% names(samples)) {
    results$reliability <- apply(samples$reliability_state, 2, mean)
  }
  
  # Reorder columns
  results <- results[, c("domain", "domain_id", "mean", "sd", 
                         "q025", "q500", "q975", "pop_share", 
                         setdiff(names(results), c("domain", "domain_id", "mean", "sd", 
                                                   "q025", "q500", "q975", "pop_share")))]
  
  # Sort by domain label
  results <- results[order(results$domain), ]
  rownames(results) <- NULL
  
  attr(results, "type") <- type
  attr(results, "prob") <- prob
  class(results) <- c("bhf_domain_estimates", "data.frame")
  
  return(results)
}


#' Extract Overall Population Estimate
#'
#' Extracts the posterior distribution of the overall (population-weighted)
#' probability from a fitted BHF model.
#'
#' @param fit An object of class \code{bhf_fit} from \code{bhf_fit()}.
#' @param prob Numeric. Probability for credible interval. Default is 0.95.
#'
#' @return A list with components:
#'   \item{mean}{Posterior mean}
#'   \item{sd}{Posterior standard deviation}
#'   \item{q025}{Lower credible interval bound}
#'   \item{q500}{Posterior median}
#'   \item{q975}{Upper credible interval bound}
#'
#' @export
overall_estimate <- function(fit, prob = 0.95) {
  
  if (!inherits(fit, "bhf_fit")) {
    stop("'fit' must be a 'bhf_fit' object from bhf_fit()", call. = FALSE)
  }
  
  # Quantile probabilities
  alpha <- (1 - prob) / 2
  probs <- c(alpha, 0.5, 1 - alpha)
  
  # Extract posterior samples
  samples <- rstan::extract(fit$stanfit)
  p_overall <- samples$p_overall
  
  list(
    mean = mean(p_overall),
    sd = sd(p_overall),
    q025 = quantile(p_overall, probs[1]),
    q500 = quantile(p_overall, probs[2]),
    q975 = quantile(p_overall, probs[3])
  )
}


#' Extract Log-Likelihood for LOO-CV
#'
#' Extracts the log-likelihood matrix for use with the loo package for
#' leave-one-out cross-validation.
#'
#' @param fit An object of class \code{bhf_fit} from \code{bhf_fit()}.
#'
#' @return A matrix of dimension (n_iterations x N) containing pointwise
#'   log-likelihood values.
#'
#' @examples
#' \dontrun{
#' # After fitting
#' fit <- bhf_fit(prepared_data, model = model)
#'
#' # Get log-likelihood and compute LOO
#' library(loo)
#' ll <- log_lik(fit)
#' loo_result <- loo(ll)
#' print(loo_result)
#' }
#'
#' @export
log_lik <- function(fit) {
  
  if (!inherits(fit, "bhf_fit")) {
    stop("'fit' must be a 'bhf_fit' object from bhf_fit()", call. = FALSE)
  }
  
  samples <- rstan::extract(fit$stanfit)
  samples$log_lik
}
