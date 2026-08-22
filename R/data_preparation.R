# Build collision-safe PSU indices for PSUs nested within strata.
#
# The construction deliberately loops over numeric stratum IDs and matches PSU
# labels only within each stratum. It never creates a pasted composite key, so
# labels containing separator-like text cannot collide.
build_nested_psu_index <- function(stratum_id, psu) {
  if (!is.numeric(stratum_id) || length(stratum_id) < 1L ||
      length(stratum_id) != length(psu) || anyNA(stratum_id) ||
      any(!is.finite(stratum_id)) ||
      any(stratum_id != as.integer(stratum_id))) {
    stop("stratum_id must be a non-empty integer vector aligned with psu",
         call. = FALSE)
  }
  if (!is.atomic(psu) || is.list(psu) || anyNA(psu)) {
    stop("psu must be an aligned atomic vector with no missing values",
         call. = FALSE)
  }

  stratum_id <- as.integer(stratum_id)
  strata <- sort(unique(stratum_id))
  if (!identical(strata, seq_len(max(stratum_id)))) {
    stop("stratum_id must contain consecutive IDs starting at 1", call. = FALSE)
  }

  psu_label <- as.character(psu)
  if (anyNA(psu_label) || any(!nzchar(trimws(psu_label)))) {
    stop("psu labels must be non-missing and non-blank", call. = FALSE)
  }

  H <- max(stratum_id)
  J_h <- integer(H)
  local <- integer(length(stratum_id))
  flat <- integer(length(stratum_id))
  mapping <- vector("list", H)
  offset <- 0L

  for (h in seq_len(H)) {
    row_index <- which(stratum_id == h)
    labels_h <- sort(unique(psu_label[row_index]), method = "radix")
    local_h <- as.integer(match(psu_label[row_index], labels_h))
    J_h[h] <- length(labels_h)
    local[row_index] <- local_h
    flat[row_index] <- offset + local_h
    mapping[[h]] <- data.frame(
      stratum_id = rep.int(h, J_h[h]),
      psu_label = labels_h,
      local = seq_len(J_h[h]),
      flat = offset + seq_len(J_h[h]),
      stringsAsFactors = FALSE
    )
    offset <- offset + J_h[h]
  }

  start <- c(1L, 1L + cumsum(J_h[-H]))

  list(
    J = as.integer(offset),
    J_h = as.integer(J_h),
    start = as.integer(start),
    local = as.integer(local),
    flat = as.integer(flat),
    mapping = do.call(rbind, mapping)
  )
}

# Scale raw survey weights for the pseudo-likelihood contract.
#
# The canonical method performs one global normalization to mean one. The
# legacy method first applies domain-specific D2 scaling and then performs the
# same final global sum-N normalization in R.
scale_likelihood_weights <- function(raw_weights,
                                     method = c("mean_one", "legacy_d2"),
                                     state_id = NULL) {
  method <- match.arg(method)
  if (!is.numeric(raw_weights) || length(raw_weights) < 1L) {
    stop("raw_weights must be a non-empty numeric vector", call. = FALSE)
  }
  if (anyNA(raw_weights) || any(!is.finite(raw_weights)) ||
      any(raw_weights <= 0)) {
    stop("raw_weights must contain only finite, strictly positive values",
         call. = FALSE)
  }

  N <- length(raw_weights)
  original_names <- names(raw_weights)

  if (identical(method, "mean_one")) {
    anchored <- raw_weights / max(raw_weights)
    scaled <- anchored / mean(anchored)
  } else {
    warning(structure(
      list(
        message = paste0(
          "weight_scaling='legacy_d2' is deprecated and retained only for ",
          "sensitivity/reproduction; removal will be no earlier than 0.5.0."
        ),
        call = NULL
      ),
      class = c("bhf_legacy_scaling_warning", "warning", "condition")
    ))
    if (is.null(state_id) || length(state_id) != N || !is.atomic(state_id) ||
        is.list(state_id) || anyNA(state_id)) {
      stop("state_id must be an aligned, non-missing atomic vector for legacy_d2",
           call. = FALSE)
    }
    state_label <- as.character(state_id)
    if (anyNA(state_label) || any(!nzchar(trimws(state_label))) ||
        (is.numeric(state_id) && any(!is.finite(state_id)))) {
      stop("state_id must contain finite, non-blank labels for legacy_d2",
           call. = FALSE)
    }

    d2_weights <- numeric(N)
    for (row_index in split(seq_len(N), state_label)) {
      weights_s <- raw_weights[row_index]
      anchored_s <- weights_s / max(weights_s)
      effective_n_s <- sum(anchored_s)^2 / sum(anchored_s^2)
      d2_weights[row_index] <- anchored_s *
        (effective_n_s / sum(anchored_s))
    }
    scaled <- d2_weights / mean(d2_weights)
  }

  names(scaled) <- original_names
  scaled
}

#' Prepare Data for BHF Model
#'
#' Transforms survey data into the format required by the BHF Stan model.
#' This function handles index recoding, weight scaling, and computation of
#' Taylor-linearized domain sampling variances needed for de-attenuation.
#'
#' @param data A data frame containing the survey data.
#' @param outcome Character string. Name of the binary outcome variable (0/1).
#' @param domain Character string. Name of the domain/state variable.
#' @param strata Character string. Name of the stratification variable.
#' @param psu Character string. Name of the PSU (primary sampling unit) variable.
#' @param weights Character string. Name of the survey weight variable.
#' @param population_shares Optional named positive numeric vector whose names
#'   exactly match the observed complete-case domains. Values are normalized to
#'   sum to one in canonical domain order. If NULL, shares are estimated from
#'   retained raw survey weights.
#' @param weight_scaling Likelihood-weight scaling method. The default
#'   \code{"mean_one"} scales raw weights globally to sum to the retained sample
#'   size. \code{"legacy_d2"} is a warned sensitivity/reproduction path.
#' @param deattenuation Character mode: `"taylor"` (default), `"supplied"`,
#'   or `"none"`.
#' @param sampling_variances Named nonnegative domain sampling variances for
#'   `deattenuation = "supplied"`.
#' @param sampling_variance_method Provenance label for supplied variances:
#'   `"external_taylor"`, `"external_replicate"`, or `"external_other"`.
#' @param use_deattenuation Deprecated logical compatibility argument. Use
#'   `deattenuation` instead.
#' @param prior_alpha_mean Numeric. Prior mean for the intercept on logit scale.
#'   Default is NULL, which estimates from data.
#' @param prior_alpha_sd Numeric. Prior SD for the intercept. Default is 0.5.
#' @param sigma_state_prior Article prior-sensitivity selector for the state
#'   random-effect SD. Exactly one of \code{"half_t3_2.5"} (baseline),
#'   \code{"half_normal_1"}, \code{"half_cauchy_2.5"}, or
#'   \code{"half_t3_5"}. Stratum and PSU SD priors remain half-t(3, 0, 2.5).
#'
#' @return An object of class \code{bhf_data} containing:
#'   \item{schema_version}{Prepared-data schema version.}
#'   \item{contract_id}{Stable prepared-data contract identifier.}
#'   \item{stan_data}{List of data formatted for Stan}
#'   \item{mapping}{List containing domain/strata/PSU label mappings}
#'   \item{row_provenance}{Original/retained/dropped row mapping and ledger.}
#'   \item{provenance}{Consolidated row, weight, share, vhat, and prior metadata.}
#'   \item{analysis_data}{Compact retained analysis frame with raw and
#'     likelihood weights.}
#'   \item{domain_summary}{Data frame with domain-level summary statistics}
#'   \item{input_info}{List recording input column names and settings}
#'
#' @details
#' This function performs several critical transformations:
#'
#' \describe{
#'   \item{Index Recoding}{All grouping variables are recoded to consecutive
#'     integers starting from 1 (required by Stan).}
#'   \item{Weight Scaling}{Positive raw weights are retained and globally
#'     normalized once in R to have mean one.}
#'   \item{Sampling Variance Estimation}{For each domain, estimates the Taylor-
#'     linearized sampling variance of the proportion from the declared
#'     stratified PSU survey design and raw weights.}
#'   \item{PSU Structure}{Creates the nested PSU-within-stratum structure
#'     required by the Stan model.}
#' }
#'
#' @section Weight Scaling:
#' The default likelihood weight is
#' \deqn{w^*_i = w_i \times \frac{N}{\sum_i w_i},}
#' so the scaled weights sum to the retained sample size N while preserving all
#' raw-weight ratios. Legacy D2 scaling is available only through the explicit
#' warned \code{weight_scaling = "legacy_d2"} sensitivity path; it receives a
#' final global sum-N normalization in R.
#'
#' @section Sampling-variance contract:
#' The default \code{deattenuation = "taylor"} uses the retained raw survey
#' weights and the one-stage stratified PSU design. Singleton strata fail
#' explicitly; no cap, floor, or default design effect is substituted.
#' \code{deattenuation = "supplied"} treats named domain variances as fixed
#' external inputs. Their estimation uncertainty is not propagated.
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
                              weight_scaling = c("mean_one", "legacy_d2"),
                              deattenuation = NULL,
                              sampling_variances = NULL,
                              sampling_variance_method = NULL,
                              use_deattenuation = NULL,
                              prior_alpha_mean = NULL,
                              prior_alpha_sd = 0.5,
                              sigma_state_prior = names(
                                sigma_state_prior_catalog()
                              )) {

  contract <- bhf_data_contract_spec()
  prior_alpha_sd_source <- if (missing(prior_alpha_sd)) {
    "default"
  } else {
    "user_supplied"
  }
  
  # ==========================================================================
  # Complete-Case Analysis Universe and Input Validation
  # ==========================================================================

  column_arguments <- list(outcome, domain, strata, psu, weights)
  valid_column_arguments <- vapply(
    column_arguments,
    function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x),
    logical(1)
  )

  # Delegate malformed top-level input to the central validator before using
  # the column specifications to construct the complete-case universe.
  if (!is.data.frame(data) || !all(valid_column_arguments)) {
    validate_input_data(data, outcome, domain, strata, psu, weights)
  }

  required_variables <- stats::setNames(
    c(outcome, domain, strata, psu, weights),
    c("outcome", "domain", "strata", "psu", "weights")
  )
  missing_variables <- setdiff(required_variables, names(data))
  if (length(missing_variables) > 0L) {
    validate_input_data(data, outcome, domain, strata, psu, weights)
  }

  n_original <- nrow(data)
  required_data <- data[, unname(required_variables), drop = FALSE]
  names(required_data) <- names(required_variables)
  missing_matrix <- is.na(required_data)
  retained_rows <- which(rowSums(missing_matrix) == 0L)
  dropped_rows <- which(rowSums(missing_matrix) > 0L)

  missing_entries <- which(missing_matrix, arr.ind = TRUE)
  if (nrow(missing_entries) > 0L) {
    missing_entries <- missing_entries[
      order(missing_entries[, "row"], missing_entries[, "col"]),
      ,
      drop = FALSE
    ]
    missing_reason_ledger <- data.frame(
      original_row = as.integer(missing_entries[, "row"]),
      original_row_name = rownames(data)[missing_entries[, "row"]],
      field = names(required_variables)[missing_entries[, "col"]],
      variable = unname(required_variables[missing_entries[, "col"]]),
      stringsAsFactors = FALSE
    )
  } else {
    missing_reason_ledger <- data.frame(
      original_row = integer(),
      original_row_name = character(),
      field = character(),
      variable = character(),
      stringsAsFactors = FALSE
    )
  }

  if (length(dropped_rows) > 0L) {
    missing_counts <- colSums(missing_matrix)
    missing_counts <- missing_counts[missing_counts > 0L]
    warning(
      "Excluded ", length(dropped_rows),
      " observations with missing required analysis values (",
      paste0(names(missing_counts), "=", missing_counts, collapse = ", "),
      ").",
      call. = FALSE
    )
  }

  data <- data[retained_rows, , drop = FALSE]
  validate_input_data(data, outcome, domain, strata, psu, weights)
  weight_scaling <- match.arg(weight_scaling)
  sigma_state_prior_info <- resolve_sigma_state_prior(sigma_state_prior)
  validate_preparation_options(
    use_deattenuation = use_deattenuation,
    prior_alpha_mean = prior_alpha_mean,
    prior_alpha_sd = prior_alpha_sd
  )
  
  # Extract columns
  df <- data.frame(
    y = data[[outcome]],
    domain = data[[domain]],
    stratum = data[[strata]],
    psu = data[[psu]],
    weight = data[[weights]],
    stringsAsFactors = FALSE
  )
  
  # ==========================================================================
  # Recode Indices (Stan requires 1-indexed consecutive integers)
  # ==========================================================================
  
  # Domain (state) mapping
  domain_labels <- as.character(sort(unique(df$domain)))
  domain_map <- stats::setNames(seq_along(domain_labels), domain_labels)
  df$state_id <- domain_map[as.character(df$domain)]
  
  # Stratum mapping
  stratum_labels <- as.character(sort(unique(df$stratum)))
  stratum_map <- stats::setNames(seq_along(stratum_labels), stratum_labels)
  df$stratum_id <- stratum_map[as.character(df$stratum)]
  
  # Dimensions
  N <- nrow(df)
  S <- length(domain_labels)
  H <- length(stratum_labels)
  
  # ==========================================================================
  # PSU Structure Within Strata
  # ==========================================================================

  psu_index <- build_nested_psu_index(df$stratum_id, df$psu)
  J <- psu_index$J
  J_h <- psu_index$J_h
  psu_start <- psu_index$start
  df$psu_in_stratum_id <- psu_index$local
  df$psu_flat_id <- psu_index$flat
  
  # ==========================================================================
  # Compute Domain-Level Statistics
  # ==========================================================================
  
  domain_summary <- compute_domain_summary(df, S)
  
  # ==========================================================================
  # Likelihood Weight Scaling
  # ==========================================================================

  df$w_lik <- scale_likelihood_weights(
    df$weight,
    method = weight_scaling,
    state_id = if (identical(weight_scaling, "legacy_d2")) df$state_id else NULL
  )

  if (identical(weight_scaling, "mean_one")) {
    global_factor <- N / sum(df$weight)
    legacy_domain_factors <- NULL
  } else {
    d2_domain_factors <- domain_summary$eff_n / domain_summary$weight_sum
    d2_intermediate <- df$weight * d2_domain_factors[df$state_id]
    global_factor <- N / sum(d2_intermediate)
    legacy_domain_factors <- stats::setNames(d2_domain_factors, domain_labels)
  }

  weight_provenance <- list(
    method = weight_scaling,
    normalization = "global_sum_N_in_R",
    N = N,
    raw_sum = sum(df$weight),
    likelihood_sum = sum(df$w_lik),
    global_factor = global_factor,
    legacy_domain_factors = legacy_domain_factors
  )
  
  # ==========================================================================
  # Population Shares
  # ==========================================================================
  
  population_share_info <- resolve_population_shares(
    population_shares = population_shares,
    domain_labels = domain_labels,
    state_id = df$state_id,
    raw_weights = df$weight
  )
  w_state_pop_share <- unname(population_share_info$values)
  
  # ==========================================================================
  # Sampling Variance for De-attenuation
  # ==========================================================================

  analysis_data <- data.frame(
    original_row = retained_rows,
    y = as.integer(df$y),
    domain_label = as.character(df$domain),
    state_id = as.integer(df$state_id),
    stratum_label = as.character(df$stratum),
    stratum_id = as.integer(df$stratum_id),
    psu_label = as.character(df$psu),
    psu_in_stratum_id = as.integer(df$psu_in_stratum_id),
    psu_flat_id = as.integer(df$psu_flat_id),
    raw_weight = as.numeric(df$weight),
    w_lik = as.numeric(df$w_lik),
    stringsAsFactors = FALSE
  )
  sampling_variance_info <- resolve_sampling_variances(
    deattenuation = deattenuation,
    sampling_variances = sampling_variances,
    sampling_variance_method = sampling_variance_method,
    domain_labels = domain_labels,
    analysis_data = analysis_data,
    use_deattenuation = use_deattenuation
  )
  vhat_state <- sampling_variance_info$stan_values
  
  # ==========================================================================
  # Prior for Intercept
  # ==========================================================================

  prior_alpha_mean_source <- if (is.null(prior_alpha_mean)) {
    "estimated_from_raw_weighted_outcome_logit"
  } else {
    "user_supplied"
  }
  design_weighted_prevalence <- sum(df$weight * df$y) / sum(df$weight)
  if (is.null(prior_alpha_mean)) {
    if (design_weighted_prevalence <= 0 || design_weighted_prevalence >= 1) {
      stop(
        "Automatic prior_alpha_mean requires a raw-weighted outcome ",
        "prevalence strictly between 0 and 1; supply prior_alpha_mean ",
        "explicitly for an all-zero or all-one outcome.",
        call. = FALSE
      )
    }
    prior_alpha_mean <- stats::qlogis(design_weighted_prevalence)
  }

  prior_info <- list(
    alpha = list(
      mean = as.numeric(prior_alpha_mean),
      mean_source = prior_alpha_mean_source,
      design_weighted_prevalence = design_weighted_prevalence,
      sd = as.numeric(prior_alpha_sd),
      sd_source = prior_alpha_sd_source,
      scale = "logit"
    ),
    random_effect_sd = list(
      family = "half_student_t",
      df = 3,
      location = 0,
      scale = 2.5,
      applies_to = c("state", "stratum", "psu")
    ),
    sigma_state = sigma_state_prior_info
  )
  
  # ==========================================================================
  # Assemble Stan Data List
  # ==========================================================================
  
  stan_data <- list(
    # Versioned data contract
    data_schema_version = contract$stan_data_schema_version,

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
    psu_flat_id = as.array(df$psu_flat_id),
    
    # Weights
    w_lik = as.array(df$w_lik),
    w_state_pop_share = as.array(w_state_pop_share),
    
    # Priors
    prior_alpha_mean = prior_alpha_mean,
    prior_alpha_sd = prior_alpha_sd,
    sigma_state_prior_code = sigma_state_prior_info$code,
    
    # De-attenuation
    use_deattenuation = as.integer(sampling_variance_info$enabled),
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
      label = psu_index$mapping$psu_label,
      id = psu_index$mapping$flat,
      stratum_label = stratum_labels[psu_index$mapping$stratum_id],
      stratum_id = psu_index$mapping$stratum_id,
      within_stratum_id = psu_index$mapping$local,
      stringsAsFactors = FALSE
    )
  )
  
  input_info <- list(
    schema_version = contract$schema_version,
    contract_id = contract$contract_id,
    stan_data_schema_version = contract$stan_data_schema_version,
    outcome = outcome,
    domain = domain,
    strata = strata,
    psu = psu,
    weights = weights,
    weight_scaling = weight_scaling,
    sigma_state_prior = sigma_state_prior_info$variant,
    population_share_source = population_share_info$source,
    deattenuation = sampling_variance_info$provenance$mode,
    use_deattenuation = sampling_variance_info$enabled,
    n_original = n_original,
    n_used = N,
    retained_rows = retained_rows,
    dropped_rows = dropped_rows,
    missing_reason_ledger = missing_reason_ledger
  )

  row_provenance <- assemble_row_provenance(
    n_original = n_original,
    n_used = N,
    retained_rows = retained_rows,
    dropped_rows = dropped_rows,
    missing_reason_ledger = missing_reason_ledger
  )

  population_share_provenance <- list(
    source = population_share_info$source,
    domain_labels = domain_labels,
    input_sum = population_share_info$input_sum,
    normalization_factor = population_share_info$normalization_factor
  )

  provenance <- assemble_data_provenance(
    contract = contract,
    rows = row_provenance,
    weights = weight_provenance,
    population_shares = population_share_provenance,
    sampling_variances = sampling_variance_info$provenance,
    prior = prior_info
  )
  
  result <- list(
    schema_version = contract$schema_version,
    contract_id = contract$contract_id,
    stan_data_schema_version = contract$stan_data_schema_version,
    stan_data = stan_data,
    mapping = mapping,
    row_provenance = row_provenance,
    provenance = provenance,
    prior_info = prior_info,
    weight_info = list(
      original_row = retained_rows,
      raw = as.numeric(df$weight),
      likelihood = as.numeric(df$w_lik),
      provenance = weight_provenance
    ),
    population_share_info = population_share_info,
    sampling_variance_info = sampling_variance_info,
    analysis_data = analysis_data,
    domain_summary = cbind(
      domain_summary,
      data.frame(
        pop_share = w_state_pop_share,
        vhat = if (sampling_variance_info$enabled) {
          vhat_state
        } else {
          rep(NA_real_, S)
        }
      )
    ),
    input_info = input_info
  )
  
  class(result) <- c("bhf_data", "list")

  validate_bhf_data_contract(result)
  
  return(result)
}


#' @export
print.bhf_data <- function(x, ...) {
  if (exists(".bhf_legacy_assert_current_bhf_data", mode = "function")) {
    .bhf_legacy_assert_current_bhf_data(x)
  } else {
    validate_bhf_data_contract(x)
  }
  cat("=== BHF Prepared Data (schema ", x$schema_version, ") ===\n\n",
      sep = "")
  cat("Rows: ", x$stan_data$N, " retained / ",
      x$row_provenance$n_original, " original",
      if (length(x$row_provenance$dropped_rows)) {
        paste0(" (", length(x$row_provenance$dropped_rows), " dropped)")
      } else "", "\n", sep = "")
  cat("Structure: ", x$stan_data$S, " domains, ", x$stan_data$H,
      " strata, ", x$stan_data$J, " nested PSUs\n", sep = "")
  cat("Likelihood weights: ", x$weight_info$provenance$method,
      "; sum = ", format(sum(x$stan_data$w_lik), digits = 8), "\n",
      sep = "")
  cat("Population shares: ", x$population_share_info$source, "\n", sep = "")
  cat("Sampling variance / A*: ",
      x$sampling_variance_info$provenance$mode,
      if (isTRUE(x$sampling_variance_info$enabled))
        " (fixed vhat; uncertainty not propagated)" else " (unavailable)",
      "\n", sep = "")
  cat("State-SD prior: ", x$prior_info$sigma_state$article_label,
      if (isTRUE(x$prior_info$sigma_state$baseline)) " [baseline]" else
        " [sensitivity]", "\n", sep = "")
  cat("Use with: fit <- bhf_fit(prepared_data, model = compiled_model)\n")
  invisible(x)
}


#' @export
summary.bhf_data <- function(object, ...) {
  if (exists(".bhf_legacy_assert_current_bhf_data", mode = "function")) {
    .bhf_legacy_assert_current_bhf_data(object)
  } else {
    validate_bhf_data_contract(object)
  }
  cat("=== BHF Prepared Data Summary ===\n\n")
  
  cat("Input Information:\n")
  cat("  Outcome variable: ", object$input_info$outcome, "\n")
  cat("  Domain variable: ", object$input_info$domain, "\n")
  cat("  Strata variable: ", object$input_info$strata, "\n")
  cat("  PSU variable: ", object$input_info$psu, "\n")
  cat("  Weight variable: ", object$input_info$weights, "\n")
  cat("  Original N: ", object$input_info$n_original, "\n")
  cat("  Used N: ", object$input_info$n_used, "\n")
  cat("  Dropped N: ", length(object$row_provenance$dropped_rows), "\n")
  
  cat("\nDimensions:\n")
  cat("  Observations: ", object$stan_data$N, "\n")
  cat("  Domains: ", object$stan_data$S, "\n")
  cat("  Strata: ", object$stan_data$H, "\n")
  cat("  PSUs: ", object$stan_data$J, "\n")

  cat("\nScientific Inputs:\n")
  cat("  Weight scaling: ", object$weight_info$provenance$method, "\n")
  cat("  Population-share source: ",
      object$population_share_info$source, "\n")
  cat("  Sampling-variance mode: ",
      object$sampling_variance_info$provenance$mode, "\n")
  cat("  Alpha prior: Normal(",
      format(object$prior_info$alpha$mean, digits = 6), ", ",
      format(object$prior_info$alpha$sd, digits = 6), ")\n", sep = "")
  cat("  Random-effect SD baseline: half-t(3, 0, 2.5); state selector: ",
      object$prior_info$sigma_state$variant, "\n", sep = "")
  
  cat("\nDomain Summary (first 10 shown):\n")
  print(utils::head(object$domain_summary, 10))
  
  invisible(object)
}
