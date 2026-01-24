#' Validate Input Data
#'
#' Internal function to validate input data for prepare_bhf_data().
#'
#' @param data Data frame to validate.
#' @param outcome Name of outcome variable.
#' @param domain Name of domain variable.
#' @param strata Name of strata variable.
#' @param psu Name of PSU variable.
#' @param weights Name of weight variable.
#'
#' @return NULL (invisibly). Throws an error if validation fails.
#'
#' @keywords internal
validate_input_data <- function(data, outcome, domain, strata, psu, weights) {
  
  # Check that data is a data frame

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }
  
  # Check required variables exist
  required_vars <- c(outcome, domain, strata, psu, weights)
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(
      "Missing required variables in data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Check outcome is binary
  y <- data[[outcome]]
  y_valid <- y[!is.na(y)]
  
  if (!all(y_valid %in% c(0, 1))) {
    stop(
      "Outcome variable '", outcome, "' must be binary (0/1).\n",
      "Found values: ", paste(unique(y_valid), collapse = ", "),
      call. = FALSE
    )
  }
  
  if (length(y_valid) == 0) {
    stop("Outcome variable '", outcome, "' has no non-missing values", call. = FALSE)
  }
  
  # Check weights are positive
  w <- data[[weights]]
  
  if (any(w <= 0, na.rm = TRUE)) {
    stop(
      "Weight variable '", weights, "' must be strictly positive.\n",
      "Found ", sum(w <= 0, na.rm = TRUE), " non-positive values.",
      call. = FALSE
    )
  }
  
  # Check sample size
  if (nrow(data) < 50) {
    warning(
      "Sample size (", nrow(data), ") is very small. ",
      "Results may be unreliable.",
      call. = FALSE
    )
  }
  
  # Check number of domains
  n_domains <- length(unique(data[[domain]]))
  if (n_domains < 3) {
    warning(
      "Only ", n_domains, " domains found. ",
      "Variance decomposition requires multiple domains.",
      call. = FALSE
    )
  }
  
  invisible(NULL)
}


#' Validate Stan Data
#'
#' Internal function to validate the Stan data list before sampling.
#'
#' @param stan_data List of Stan data.
#'
#' @return NULL (invisibly). Throws an error if validation fails.
#'
#' @keywords internal
validate_stan_data <- function(stan_data) {
  
  # Check dimensions
  stopifnot(stan_data$N > 0)
  stopifnot(stan_data$S > 0)
  stopifnot(stan_data$H > 0)
  stopifnot(stan_data$J > 0)
  
  # Check index ranges
  if (any(stan_data$state_id < 1) || any(stan_data$state_id > stan_data$S)) {
    stop("state_id values must be in range [1, S]", call. = FALSE)
  }
  
  if (any(stan_data$stratum_id < 1) || any(stan_data$stratum_id > stan_data$H)) {
    stop("stratum_id values must be in range [1, H]", call. = FALSE)
  }
  
  # Check lengths
  if (length(stan_data$y) != stan_data$N) {
    stop("Length of y must equal N", call. = FALSE)
  }
  
  if (length(stan_data$state_id) != stan_data$N) {
    stop("Length of state_id must equal N", call. = FALSE)
  }
  
  if (length(stan_data$w_lik) != stan_data$N) {
    stop("Length of w_lik must equal N", call. = FALSE)
  }
  
  if (length(stan_data$J_h) != stan_data$H) {
    stop("Length of J_h must equal H", call. = FALSE)
  }
  
  if (sum(stan_data$J_h) != stan_data$J) {
    stop("Sum of J_h must equal J", call. = FALSE)
  }
  
  # Check population shares sum to 1
  if (abs(sum(stan_data$w_state_pop_share) - 1) > 1e-6) {
    stop("w_state_pop_share must sum to 1", call. = FALSE)
  }
  
  # Check for NA/NaN/Inf
  numeric_fields <- c("y", "w_lik", "w_state_pop_share", "vhat_state")
  for (field in numeric_fields) {
    vals <- stan_data[[field]]
    if (any(is.na(vals)) || any(is.nan(vals)) || any(is.infinite(vals))) {
      stop("Stan data field '", field, "' contains NA, NaN, or Inf values", call. = FALSE)
    }
  }
  
  invisible(NULL)
}


#' Compute Domain Summary Statistics
#'
#' Internal function to compute summary statistics for each domain.
#'
#' @param df Data frame with columns: y, state_id, weight.
#' @param S Number of domains.
#'
#' @return Data frame with domain-level statistics.
#'
#' @keywords internal
compute_domain_summary <- function(df, S) {
  
  # Initialize
  summary_list <- vector("list", S)
  
  for (s in seq_len(S)) {
    idx <- which(df$state_id == s)
    
    w_s <- df$weight[idx]
    y_s <- df$y[idx]
    
    n_s <- length(idx)
    n_obs_s <- sum(!is.na(y_s))
    
    # Weighted proportion
    y_valid <- y_s[!is.na(y_s)]
    w_valid <- w_s[!is.na(y_s)]
    
    if (n_obs_s > 0) {
      prop_weighted <- sum(w_valid * y_valid) / sum(w_valid)
    } else {
      prop_weighted <- NA_real_
    }
    
    # Effective sample size
    eff_n <- (sum(w_s)^2) / sum(w_s^2)
    
    summary_list[[s]] <- data.frame(
      state_id = s,
      n = n_s,
      n_obs = n_obs_s,
      prop_weighted = prop_weighted,
      weight_sum = sum(w_s),
      eff_n = eff_n,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, summary_list)
}


#' Scale Weights Using Method D2
#'
#' Internal function to scale survey weights using the D2 method
#' (Pfeffermann et al., 1998).
#'
#' @param weight Original survey weights.
#' @param state_id State/domain IDs.
#' @param domain_summary Domain summary from compute_domain_summary().
#'
#' @return Scaled weights.
#'
#' @keywords internal
scale_weights_d2 <- function(weight, state_id, domain_summary) {
  
  w_scaled <- numeric(length(weight))
  
  for (s in seq_len(nrow(domain_summary))) {
    idx <- which(state_id == s)
    w_s <- weight[idx]
    
    eff_n_s <- domain_summary$eff_n[s]
    weight_sum_s <- domain_summary$weight_sum[s]
    
    # Method D2: scale so sum of scaled weights equals effective n
    w_scaled[idx] <- w_s * (eff_n_s / weight_sum_s)
  }
  
  return(w_scaled)
}


#' Estimate Sampling Variance for Each Domain
#'
#' Internal function to estimate the sampling variance of the
#' domain-specific proportion.
#'
#' @param domain_summary Domain summary from compute_domain_summary().
#'
#' @return Numeric vector of sampling variance estimates.
#'
#' @keywords internal
estimate_sampling_variance <- function(domain_summary) {
  
  S <- nrow(domain_summary)
  vhat <- numeric(S)
  
  for (s in seq_len(S)) {
    p_s <- domain_summary$prop_weighted[s]
    n_s <- domain_summary$n_obs[s]
    eff_n_s <- domain_summary$eff_n[s]
    
    # Skip if no observations
    if (is.na(p_s) || n_s == 0) {
      vhat[s] <- 0
      next
    }
    
    # Design effect estimate
    if (eff_n_s > 0 && n_s > 0) {
      deff_s <- n_s / eff_n_s
    } else {
      deff_s <- 1.5  # Conservative default
    }
    
    # Bound design effect
    deff_s <- pmax(1, pmin(10, deff_s))
    
    # Sampling variance estimate
    # V(p) = deff * p(1-p) / n
    p_s_bounded <- pmax(0.01, pmin(0.99, p_s))
    
    if (n_s > 1) {
      vhat[s] <- deff_s * p_s_bounded * (1 - p_s_bounded) / n_s
    } else {
      vhat[s] <- 0.1  # Conservative default for small samples
    }
  }
  
  return(vhat)
}


#' Calculate Effective Sample Size
#'
#' Calculates the effective sample size given survey weights.
#'
#' @param weights Numeric vector of survey weights.
#'
#' @return Effective sample size (scalar).
#'
#' @examples
#' w <- c(1, 2, 1.5, 3, 2)
#' calc_eff_n(w)  # Should be less than length(w)
#'
#' @export
calc_eff_n <- function(weights) {
  w <- as.numeric(weights)
  if (any(w < 0)) {
    stop("Weights must be non-negative", call. = FALSE)
  }
  (sum(w)^2) / sum(w^2)
}
