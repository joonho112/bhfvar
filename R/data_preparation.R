#' Prepare Data for BHF Model
#'
#' Transforms survey data into the format required by the BHF Stan model.
#' This function handles index recoding, weight scaling, and computation
#' of design effect estimates needed for de-attenuation.
#'
#' @param data A data frame containing the survey data.
#' @param outcome Character string. Name of the binary outcome variable (0/1).
#' @param domain Character string. Name of the domain/state variable.
#' @param strata Character string. Name of the stratification variable.
#' @param psu Character string. Name of the PSU (primary sampling unit) variable.
#' @param weights Character string. Name of the survey weight variable.
#' @param population_shares Optional numeric vector of length S (number of domains)
#'   containing population shares for each domain. Must sum to 1. If NULL,
#'   shares are estimated from the weighted data.
#' @param use_deattenuation Logical. If TRUE (default), computes and applies
#'   de-attenuation adjustment for finite-sample variance inflation.
#' @param prior_alpha_mean Numeric. Prior mean for the intercept on logit scale.
#'   Default is NULL, which estimates from data.
#' @param prior_alpha_sd Numeric. Prior SD for the intercept. Default is 1.5.
#'
#' @return An object of class \code{bhf_data} containing:
#'   \item{stan_data}{List of data formatted for Stan}
#'   \item{mapping}{List containing domain/strata/PSU label mappings}
#'   \item{domain_summary}{Data frame with domain-level summary statistics}
#'   \item{input_info}{List recording input column names and settings}
#'
#' @details
#' This function performs several critical transformations:
#'
#' \describe{
#'   \item{Index Recoding}{All grouping variables are recoded to consecutive
#'     integers starting from 1 (required by Stan).}
#'   \item{Weight Scaling}{Weights are scaled using Method D2 from Pfeffermann
#'     et al. (1998) to have effective sample size within each domain.
#'     This is critical for proper pseudo-likelihood estimation.}
#'   \item{Sampling Variance Estimation}{For each domain, estimates the
#'     sampling variance of the proportion using the design effect.}
#'   \item{PSU Structure}{Creates the nested PSU-within-stratum structure
#'     required by the Stan model.}
#' }
#'
#' @section Weight Scaling (Method D2):
#' Weights are scaled so that for each domain s:
#' \deqn{w^*_i = w_i \times \frac{n^{eff}_s}{\sum_{i \in s} w_i}}
#' where \eqn{n^{eff}_s = (\sum w_i)^2 / \sum w_i^2} is the effective sample size.
#' This ensures the pseudo-likelihood contributes appropriate information.
#'
#' @section Sampling Variance Estimation:
#' The estimated sampling variance for domain s is:
#' \deqn{\hat{V}_s = \frac{deff_s \times \hat{p}_s(1-\hat{p}_s)}{n_s}}
#' where \eqn{deff_s} is the design effect. A default value of 1.5 is used
#' when the design effect cannot be reliably estimated.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(bhf_synthetic_data)
#'
#' # Prepare data for Stan
#' prepared <- prepare_bhf_data(
#'   data = bhf_synthetic_data,
#'   outcome = "has_subsidy",
#'   domain = "state",
#'   strata = "stratum",
#'   psu = "psu",
#'   weights = "weight"
#' )
#'
#' # Inspect the result
#' print(prepared)
#' summary(prepared)
#' }
#'
#' @export
prepare_bhf_data <- function(data,
                              outcome,
                              domain,
                              strata,
                              psu,
                              weights,
                              population_shares = NULL,
                              use_deattenuation = TRUE,
                              prior_alpha_mean = NULL,
                              prior_alpha_sd = 1.5) {
  
  # ==========================================================================
  # Input Validation
  # ==========================================================================
  
  validate_input_data(data, outcome, domain, strata, psu, weights)
  
  # Extract columns
  df <- data.frame(
    y = data[[outcome]],
    domain = data[[domain]],
    stratum = data[[strata]],
    psu = data[[psu]],
    weight = data[[weights]],
    stringsAsFactors = FALSE
  )
  
  # Check for missing values in design variables
  design_complete <- complete.cases(df[, c("domain", "stratum", "psu", "weight")])
  if (any(!design_complete)) {
    n_missing <- sum(!design_complete)
    warning(
      n_missing, " observations have missing values in design variables and will be excluded.",
      call. = FALSE
    )
    df <- df[design_complete, ]
  }
  
  # ==========================================================================
  # Recode Indices (Stan requires 1-indexed consecutive integers)
  # ==========================================================================
  
  # Domain (state) mapping
  domain_labels <- sort(unique(df$domain))
  domain_map <- setNames(seq_along(domain_labels), domain_labels)
  df$state_id <- domain_map[as.character(df$domain)]
  
  # Stratum mapping
  stratum_labels <- sort(unique(df$stratum))
  stratum_map <- setNames(seq_along(stratum_labels), stratum_labels)
  df$stratum_id <- stratum_map[as.character(df$stratum)]
  
  # PSU mapping (globally unique)
  psu_labels <- sort(unique(df$psu))
  psu_map <- setNames(seq_along(psu_labels), psu_labels)
  df$psu_global_id <- psu_map[as.character(df$psu)]
  
  # Dimensions
  N <- nrow(df)
  S <- length(domain_labels)
  H <- length(stratum_labels)
  J <- length(psu_labels)
  
  # ==========================================================================
  # PSU Structure Within Strata
  # ==========================================================================
  
  # Get PSUs within each stratum
  psu_stratum <- unique(df[, c("stratum_id", "psu_global_id")])
  psu_stratum <- psu_stratum[order(psu_stratum$stratum_id, psu_stratum$psu_global_id), ]
  
  # J_h: number of PSUs in each stratum
  J_h <- as.integer(table(psu_stratum$stratum_id))
  
  # psu_start: starting index for PSUs in each stratum
  psu_start <- c(1L, cumsum(J_h[-H]) + 1L)
  
  # Create PSU-within-stratum index
  psu_within_stratum <- integer(J)
  for (h in seq_len(H)) {
    idx <- which(psu_stratum$stratum_id == h)
    psu_within_stratum[psu_stratum$psu_global_id[idx]] <- seq_along(idx)
  }
  df$psu_in_stratum_id <- psu_within_stratum[df$psu_global_id]
  
  # ==========================================================================
  # Compute Domain-Level Statistics
  # ==========================================================================
  
  domain_summary <- compute_domain_summary(df, S)
  
  # ==========================================================================
  # Weight Scaling (Method D2)
  # ==========================================================================
  
  df$w_lik <- scale_weights_d2(df$weight, df$state_id, domain_summary)
  
  # ==========================================================================
  # Population Shares
  # ==========================================================================
  
  if (is.null(population_shares)) {
    # Estimate from weighted data
    w_state_pop_share <- domain_summary$weight_sum / sum(domain_summary$weight_sum)
  } else {
    if (length(population_shares) != S) {
      stop(
        "population_shares must have length ", S, " (number of domains), ",
        "but has length ", length(population_shares),
        call. = FALSE
      )
    }
    if (abs(sum(population_shares) - 1) > 1e-6) {
      stop("population_shares must sum to 1", call. = FALSE)
    }
    w_state_pop_share <- population_shares
  }
  
  # ==========================================================================
  # Sampling Variance Estimation for De-attenuation
  # ==========================================================================
  
  if (use_deattenuation) {
    vhat_state <- estimate_sampling_variance(domain_summary)
  } else {
    vhat_state <- rep(0, S)
  }
  
  # ==========================================================================
  # Prior for Intercept
  # ==========================================================================
  
  if (is.null(prior_alpha_mean)) {
    # Estimate from weighted overall proportion
    overall_prop <- sum(df$weight * df$y, na.rm = TRUE) / sum(df$weight[!is.na(df$y)])
    overall_prop <- pmax(0.01, pmin(0.99, overall_prop))  # Bound away from 0/1
    prior_alpha_mean <- qlogis(overall_prop)  # logit transform
  }
  
  # ==========================================================================
  # Assemble Stan Data List
  # ==========================================================================
  
  stan_data <- list(
    # Dimensions
    N = N,
    S = S,
    H = H,
    J = J,
    
    # Outcome
    y = as.array(as.integer(df$y)),
    
    # Indices
    state_id = as.array(df$state_id),
    stratum_id = as.array(df$stratum_id),
    
    # PSU structure
    J_h = as.array(J_h),
    psu_start = as.array(psu_start),
    psu_in_stratum_id = as.array(df$psu_in_stratum_id),
    
    # Weights
    w_lik = as.array(df$w_lik),
    w_state_pop_share = as.array(w_state_pop_share),
    
    # Priors
    prior_alpha_mean = prior_alpha_mean,
    prior_alpha_sd = prior_alpha_sd,
    
    # De-attenuation
    use_deattenuation = as.integer(use_deattenuation),
    vhat_state = as.array(vhat_state)
  )
  
  # ==========================================================================
  # Validate Stan Data
  # ==========================================================================
  
  validate_stan_data(stan_data)
  
  # ==========================================================================
  # Create Return Object
  # ==========================================================================
  
  mapping <- list(
    domain = data.frame(
      label = domain_labels,
      id = seq_along(domain_labels),
      stringsAsFactors = FALSE
    ),
    stratum = data.frame(
      label = stratum_labels,
      id = seq_along(stratum_labels),
      stringsAsFactors = FALSE
    ),
    psu = data.frame(
      label = psu_labels,
      id = seq_along(psu_labels),
      stringsAsFactors = FALSE
    )
  )
  
  input_info <- list(
    outcome = outcome,
    domain = domain,
    strata = strata,
    psu = psu,
    weights = weights,
    use_deattenuation = use_deattenuation,
    n_original = nrow(data),
    n_used = N
  )
  
  result <- list(
    stan_data = stan_data,
    mapping = mapping,
    domain_summary = cbind(
      domain_summary,
      data.frame(
        pop_share = w_state_pop_share,
        vhat = vhat_state
      )
    ),
    input_info = input_info
  )
  
  class(result) <- c("bhf_data", "list")
  
  return(result)
}


#' @export
print.bhf_data <- function(x, ...) {
  cat("=== BHF Prepared Data ===\n\n")
  cat("Sample size: ", x$stan_data$N, "\n")
  cat("Number of domains: ", x$stan_data$S, "\n")
  cat("Number of strata: ", x$stan_data$H, "\n")
  cat("Number of PSUs: ", x$stan_data$J, "\n")
  cat("\nOutcome proportion: ", 
      sprintf("%.1f%%", mean(x$stan_data$y) * 100), "\n")
  cat("De-attenuation: ", 
      ifelse(x$stan_data$use_deattenuation == 1, "enabled", "disabled"), "\n")
  cat("\nUse with: fit <- bhf_fit(prepared_data, model = compiled_model)\n")
  invisible(x)
}


#' @export
summary.bhf_data <- function(object, ...) {
  cat("=== BHF Prepared Data Summary ===\n\n")
  
  cat("Input Information:\n")
  cat("  Outcome variable: ", object$input_info$outcome, "\n")
  cat("  Domain variable: ", object$input_info$domain, "\n")
  cat("  Strata variable: ", object$input_info$strata, "\n")
  cat("  PSU variable: ", object$input_info$psu, "\n")
  cat("  Weight variable: ", object$input_info$weights, "\n")
  cat("  Original N: ", object$input_info$n_original, "\n")
  cat("  Used N: ", object$input_info$n_used, "\n")
  
  cat("\nDimensions:\n")
  cat("  Observations: ", object$stan_data$N, "\n")
  cat("  Domains: ", object$stan_data$S, "\n")
  cat("  Strata: ", object$stan_data$H, "\n")
  cat("  PSUs: ", object$stan_data$J, "\n")
  
  cat("\nDomain Summary (first 10 shown):\n")
  print(head(object$domain_summary, 10))
  
  invisible(object)
}
