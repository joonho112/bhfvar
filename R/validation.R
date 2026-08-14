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
  validate_input_spec(data, outcome, domain, strata, psu, weights)
  validate_analysis_data(data, outcome, domain, strata, psu, weights)
  invisible(NULL)
}

# Validate only the data-frame and column-role specification.  This helper is
# deliberately separate so prepare_bhf_data() can apply its one complete-case
# filter before validating retained analysis values.
validate_input_spec <- function(data, outcome, domain, strata, psu, weights) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }

  roles <- list(
    outcome = outcome,
    domain = domain,
    strata = strata,
    psu = psu,
    weights = weights
  )

  invalid <- vapply(roles, function(x) {
    !is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)
  }, logical(1))
  if (any(invalid)) {
    stop(
      "Column arguments must each be one non-empty character string: ",
      paste(names(roles)[invalid], collapse = ", "),
      call. = FALSE
    )
  }

  required_vars <- unname(unlist(roles, use.names = FALSE))
  if (anyDuplicated(required_vars)) {
    stop("outcome, domain, strata, psu, and weights must name distinct columns",
         call. = FALSE)
  }

  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0L) {
    stop(
      "Missing required variables in data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(NULL)
}

# Validate values after the complete-case analysis universe has been selected.
validate_analysis_data <- function(data, outcome, domain, strata, psu, weights) {
  if (nrow(data) < 1L) {
    stop("Analysis data must contain at least one retained row", call. = FALSE)
  }

  y <- data[[outcome]]
  if (!(is.numeric(y) || is.integer(y) || is.logical(y))) {
    stop("Outcome variable '", outcome, "' must be numeric or logical binary (0/1)",
         call. = FALSE)
  }
  if (anyNA(y) || any(!is.finite(as.numeric(y)))) {
    stop("Outcome variable '", outcome, "' must be finite and non-missing",
         call. = FALSE)
  }
  if (!all(y %in% c(0, 1))) {
    stop(
      "Outcome variable '", outcome, "' must be binary (0/1).\n",
      "Found values: ", paste(unique(y), collapse = ", "),
      call. = FALSE
    )
  }

  w <- data[[weights]]
  if (!is.numeric(w)) {
    stop("Weight variable '", weights, "' must be numeric", call. = FALSE)
  }
  if (anyNA(w) || any(!is.finite(w))) {
    stop("Weight variable '", weights, "' must be finite and non-missing",
         call. = FALSE)
  }
  if (any(w <= 0)) {
    stop(
      "Weight variable '", weights, "' must be strictly positive.\n",
      "Found ", sum(w <= 0), " non-positive values.",
      call. = FALSE
    )
  }

  for (column in c(domain, strata, psu)) {
    values <- data[[column]]
    if (!is.atomic(values) || is.list(values)) {
      stop("Grouping variable '", column, "' must be an atomic vector",
           call. = FALSE)
    }
    labels <- as.character(values)
    if (anyNA(values) || anyNA(labels) || any(!nzchar(trimws(labels)))) {
      stop("Grouping variable '", column,
           "' must have non-missing, non-blank labels", call. = FALSE)
    }
  }

  if (nrow(data) < 50L) {
    warning(
      "Sample size (", nrow(data), ") is very small. ",
      "Results may be unreliable.",
      call. = FALSE
    )
  }

  n_domains <- length(unique(data[[domain]]))
  if (n_domains < 3L) {
    warning(
      "Only ", n_domains, " domains found. ",
      "Variance decomposition requires multiple domains.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

validate_preparation_options <- function(use_deattenuation,
                                         prior_alpha_mean,
                                         prior_alpha_sd) {
  if (!is.null(use_deattenuation) &&
      (!is.logical(use_deattenuation) || length(use_deattenuation) != 1L ||
       is.na(use_deattenuation))) {
    stop("use_deattenuation must be one non-missing logical value", call. = FALSE)
  }
  if (!is.null(prior_alpha_mean) &&
      (!is.numeric(prior_alpha_mean) || length(prior_alpha_mean) != 1L ||
       is.na(prior_alpha_mean) || !is.finite(prior_alpha_mean))) {
    stop("prior_alpha_mean must be NULL or one finite numeric value", call. = FALSE)
  }
  if (!is.numeric(prior_alpha_sd) || length(prior_alpha_sd) != 1L ||
      is.na(prior_alpha_sd) || !is.finite(prior_alpha_sd) ||
      prior_alpha_sd <= 0) {
    stop("prior_alpha_sd must be one finite, strictly positive numeric value",
         call. = FALSE)
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
  required <- c(
    "data_schema_version", "N", "S", "H", "J", "y", "state_id", "stratum_id",
    "J_h", "psu_start", "psu_in_stratum_id", "psu_flat_id", "w_lik",
    "w_state_pop_share", "prior_alpha_mean", "prior_alpha_sd",
    "sigma_state_prior_code",
    "use_deattenuation", "vhat_state"
  )
  missing_fields <- setdiff(required, names(stan_data))
  if (length(missing_fields)) {
    stop("Stan data is missing required fields: ",
         paste(missing_fields, collapse = ", "), call. = FALSE)
  }

  schema_version <- stan_data$data_schema_version
  if (!is.numeric(schema_version) || length(schema_version) != 1L ||
      is.na(schema_version) || !is.finite(schema_version) ||
      schema_version != 1L) {
    stop("data_schema_version must be the integer 1", call. = FALSE)
  }

  for (field in c("N", "S", "H", "J")) {
    value <- stan_data[[field]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 1 || value != as.integer(value)) {
      stop("Stan dimension '", field, "' must be one positive integer",
           call. = FALSE)
    }
  }

  N <- as.integer(stan_data$N)
  S <- as.integer(stan_data$S)
  H <- as.integer(stan_data$H)
  J <- as.integer(stan_data$J)

  code <- stan_data$sigma_state_prior_code
  if (!is.numeric(code) || length(code) != 1L || is.na(code) ||
      !is.finite(code) || code != as.integer(code) ||
      !code %in% 1:4) {
    stop("sigma_state_prior_code must be one integer in [1, 4]",
         call. = FALSE)
  }

  expected_lengths <- c(
    y = N, state_id = N, stratum_id = N, psu_in_stratum_id = N,
    psu_flat_id = N,
    w_lik = N, J_h = H, psu_start = H,
    w_state_pop_share = S, vhat_state = S
  )
  for (field in names(expected_lengths)) {
    if (length(stan_data[[field]]) != expected_lengths[[field]]) {
      stop("Length of ", field, " must equal ", expected_lengths[[field]],
           call. = FALSE)
    }
  }

  finite_fields <- c(
    "y", "state_id", "stratum_id", "psu_in_stratum_id", "psu_flat_id", "J_h",
    "psu_start", "w_lik", "w_state_pop_share", "vhat_state"
  )
  for (field in finite_fields) {
    values <- stan_data[[field]]
    if (!is.numeric(values) || anyNA(values) || any(!is.finite(values))) {
      stop("Stan data field '", field, "' must be numeric, finite, and non-missing",
           call. = FALSE)
    }
  }

  check_integer_ids <- function(values, field, lower, upper) {
    if (any(values != as.integer(values)) || any(values < lower) ||
        any(values > upper)) {
      stop(field, " values must be integers in range [", lower, ", ", upper,
           "]", call. = FALSE)
    }
  }
  check_integer_ids(stan_data$state_id, "state_id", 1L, S)
  check_integer_ids(stan_data$stratum_id, "stratum_id", 1L, H)

  if (any(stan_data$y != as.integer(stan_data$y)) ||
      any(!stan_data$y %in% c(0, 1))) {
    stop("y values must be binary integers (0/1)", call. = FALSE)
  }
  if (any(stan_data$w_lik <= 0)) {
    stop("w_lik values must be strictly positive", call. = FALSE)
  }
  if (abs(sum(stan_data$w_lik) - N) > 1e-10 * max(1, N)) {
    stop("w_lik must use the canonical global sum-N normalization",
         call. = FALSE)
  }

  if (any(stan_data$J_h != as.integer(stan_data$J_h)) ||
      any(stan_data$J_h < 1L) || sum(stan_data$J_h) != J) {
    stop("J_h must contain positive integers summing to J", call. = FALSE)
  }
  cumulative_j <- cumsum(as.integer(stan_data$J_h))
  expected_start <- c(1L, cumulative_j[seq_len(max(0L, H - 1L))] + 1L)
  if (!identical(as.integer(stan_data$psu_start), expected_start)) {
    stop("psu_start must identify consecutive stratum PSU blocks", call. = FALSE)
  }
  local_max <- stan_data$J_h[as.integer(stan_data$stratum_id)]
  if (any(stan_data$psu_in_stratum_id !=
          as.integer(stan_data$psu_in_stratum_id)) ||
      any(stan_data$psu_in_stratum_id < 1L) ||
      any(stan_data$psu_in_stratum_id > local_max)) {
    stop("psu_in_stratum_id values must be valid within their strata",
         call. = FALSE)
  }
  check_integer_ids(stan_data$psu_flat_id, "psu_flat_id", 1L, J)
  reconstructed_flat <- stan_data$psu_start[as.integer(stan_data$stratum_id)] +
    as.integer(stan_data$psu_in_stratum_id) - 1L
  if (!identical(as.integer(stan_data$psu_flat_id),
                 as.integer(reconstructed_flat))) {
    stop("psu_flat_id must equal its stratum-block reconstruction",
         call. = FALSE)
  }

  shares <- stan_data$w_state_pop_share
  if (any(shares <= 0) || abs(sum(shares) - 1) > 1e-12) {
    stop("w_state_pop_share must be strictly positive and sum to 1",
         call. = FALSE)
  }
  if (any(stan_data$vhat_state < 0)) {
    stop("vhat_state values must be nonnegative", call. = FALSE)
  }

  validate_preparation_options(
    use_deattenuation = as.logical(stan_data$use_deattenuation),
    prior_alpha_mean = stan_data$prior_alpha_mean,
    prior_alpha_sd = stan_data$prior_alpha_sd
  )
  if (length(stan_data$use_deattenuation) != 1L ||
      !stan_data$use_deattenuation %in% c(0L, 1L)) {
    stop("use_deattenuation must be the integer 0 or 1", call. = FALSE)
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
  if (!is.numeric(weights) || length(weights) < 1L) {
    stop("weights must be a non-empty numeric vector", call. = FALSE)
  }
  if (anyNA(weights) || any(!is.finite(weights))) {
    stop("weights must contain only finite, non-missing values", call. = FALSE)
  }
  if (any(weights <= 0)) {
    stop("weights must be strictly positive", call. = FALSE)
  }
  result <- (sum(weights)^2) / sum(weights^2)
  if (!is.finite(result) || result <= 0) {
    stop("effective sample size must be finite and positive", call. = FALSE)
  }
  result
}
