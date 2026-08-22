# Shared result contract for the article-aligned 0.5.0 extractors.

.bhf_result_abort <- function(class, message, ...) {
  stop(structure(
    c(list(message = message, call = NULL), list(...)),
    class = c(class, "error", "condition")
  ))
}

.bhf_result_warn <- function(class, message, ...) {
  warning(structure(
    c(list(message = message, call = NULL), list(...)),
    class = c(class, "warning", "condition")
  ))
}

.bhf_validate_result_fit <- function(fit) {
  if (!inherits(fit, "bhf_fit")) {
    .bhf_result_abort(
      "bhf_fit_class_error",
      "'fit' must be a 'bhf_fit' object from bhf_fit()."
    )
  }
  if (exists("detect_bhf_object_schema", mode = "function") &&
      exists(".bhf_legacy_abort", mode = "function")) {
    detection <- detect_bhf_object_schema(fit)
    if (!identical(detection$object_type, "bhf_fit") ||
        !identical(detection$status, "current")) {
      .bhf_legacy_abort(detection, "bhf_fit")
    }
  }
  fit_contract_ok <- identical(fit$schema_version, "0.5.0") &&
    identical(fit$contract_id, "bhfvar-fit-contract-0.5.0")
  data_contract_ok <- is.list(fit$data) &&
    identical(fit$data$schema_version, "0.5.0") &&
    identical(fit$data$contract_id, "bhfvar-data-contract-0.5.0") &&
    identical(fit$data$stan_data_schema_version, 1L) &&
    is.list(fit$data$stan_data) &&
    identical(fit$data$stan_data$data_schema_version, 1L)
  if (!fit_contract_ok || !data_contract_ok) {
    .bhf_result_abort(
      "bhf_fit_schema_error",
      paste0(
        "Article-aligned extractors require the exact bhf_fit/data 0.5.0 ",
        "contract markers; legacy, unversioned, or masquerading fits must ",
        "be re-prepared and refit."
      ),
      observed_schema = fit$schema_version,
      observed_contract = fit$contract_id
    )
  }
  if (is.null(fit$stanfit)) {
    .bhf_result_abort(
      "bhf_fit_schema_error",
      "The bhf_fit object is missing its Stan fit or prepared-data contract."
    )
  }
  invisible(fit)
}

.bhf_match_choice <- function(value, choices, name) {
  if (missing(value) || identical(value, choices)) value <- choices[[1L]]
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !value %in% choices) {
    .bhf_result_abort(
      "bhf_argument_error",
      paste0("'", name, "' must be exactly one of: ",
             paste(shQuote(choices), collapse = ", "), "."),
      argument = name
    )
  }
  value
}

.bhf_percent_label <- function(x) {
  paste0(trimws(formatC(100 * x, format = "fg", digits = 12,
                        drop0trailing = TRUE)), "%")
}

.bhf_interval_contract <- function(prob = 0.95) {
  if (!is.numeric(prob) || length(prob) != 1L || is.na(prob) ||
      !is.finite(prob) || prob <= 0 || prob >= 1) {
    .bhf_result_abort(
      "bhf_interval_error",
      "'prob' must be one finite numeric value strictly between 0 and 1.",
      argument = "prob"
    )
  }
  tail <- (1 - prob) / 2
  quantiles <- c(lower = tail, median = 0.5, upper = 1 - tail)
  labels <- stats::setNames(vapply(quantiles, .bhf_percent_label,
                                   character(1)), names(quantiles))
  structure(
    list(
      prob = as.numeric(prob),
      quantiles = quantiles,
      labels = labels,
      interval_label = paste(.bhf_percent_label(prob), "credible interval")
    ),
    class = c("bhf_interval_contract", "list")
  )
}

.bhf_draw_summary <- function(x, interval, name = "draws") {
  if (!is.numeric(x) || !length(x) || anyNA(x) || any(!is.finite(x))) {
    .bhf_result_abort(
      "bhf_draw_contract_error",
      paste0("Posterior field '", name,
             "' must contain non-empty finite numeric draws."),
      field = name
    )
  }
  x <- as.numeric(x)
  quantiles <- stats::quantile(
    x, probs = unname(interval$quantiles), names = FALSE, type = 7
  )
  c(
    mean = mean(x),
    sd = if (length(x) > 1L) stats::sd(x) else 0,
    lower = quantiles[[1L]],
    median = quantiles[[2L]],
    upper = quantiles[[3L]],
    prob = interval$prob
  )
}

.bhf_draw_vector <- function(x, name) {
  dimensions <- dim(x)
  if (!is.null(dimensions) &&
      (length(dimensions) > 2L ||
       (length(dimensions) == 2L && dimensions[[2L]] != 1L))) {
    .bhf_result_abort(
      "bhf_draw_shape_error",
      paste0("Posterior field '", name, "' must be one scalar per draw."),
      field = name
    )
  }
  if (!is.numeric(x) || !length(x) || anyNA(x) || any(!is.finite(x))) {
    .bhf_result_abort(
      "bhf_draw_contract_error",
      paste0("Posterior field '", name,
             "' must contain non-empty finite numeric draws."),
      field = name
    )
  }
  as.numeric(x)
}

.bhf_draw_matrix <- function(x, n_columns, name) {
  if (is.null(dim(x)) && identical(n_columns, 1L)) {
    x <- matrix(x, ncol = 1L)
  }
  if (!is.matrix(x) || !is.numeric(x) || ncol(x) != n_columns ||
      nrow(x) < 1L || anyNA(x) || any(!is.finite(x))) {
    .bhf_result_abort(
      "bhf_draw_shape_error",
      paste0("Posterior field '", name, "' must be a finite draw-by-",
             n_columns, " matrix."),
      field = name
    )
  }
  x
}

.bhf_close <- function(x, y) {
  length(x) == length(y) && all(
    is.finite(x) & is.finite(y) &
      abs(x - y) <= 1e-10 + 1e-8 * pmax(abs(x), abs(y))
  )
}

# Kept as a single wrapper so focused tests and alternative fitted-object
# backends can replace only the extraction boundary.
.bhf_extract_draws <- function(fit, pars) {
  rstan::extract(fit$stanfit, pars = pars, permuted = TRUE)
}

.bhf_extract_required <- function(fit, pars) {
  draws <- .bhf_extract_draws(fit, pars)
  if (!is.list(draws)) {
    .bhf_result_abort(
      "bhf_draw_contract_error",
      "The posterior extraction backend must return a named list."
    )
  }
  missing_fields <- setdiff(pars, names(draws))
  if (length(missing_fields)) {
    .bhf_result_abort(
      "bhf_draw_contract_error",
      paste0("The fit is missing required posterior fields: ",
             paste(missing_fields, collapse = ", "), "."),
      missing_fields = missing_fields
    )
  }
  draws[pars]
}

.bhf_component_table <- function(draws, mapping, interval, section,
                                 estimand, scale, available = TRUE) {
  rows <- lapply(names(mapping), function(component) {
    field <- unname(mapping[[component]])
    values <- .bhf_draw_vector(draws[[field]], field)
    summary <- .bhf_draw_summary(values, interval, field)
    data.frame(
      section = section,
      estimand = estimand,
      scale = scale,
      component = component,
      mean = unname(summary[["mean"]]),
      sd = unname(summary[["sd"]]),
      lower = unname(summary[["lower"]]),
      median = unname(summary[["median"]]),
      upper = unname(summary[["upper"]]),
      prob = interval$prob,
      available = available,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.bhf_unavailable_table <- function(components, interval, section, estimand,
                                   scale) {
  data.frame(
    section = section,
    estimand = estimand,
    scale = scale,
    component = components,
    mean = NA_real_, sd = NA_real_, lower = NA_real_, median = NA_real_,
    upper = NA_real_, prob = interval$prob, available = FALSE,
    stringsAsFactors = FALSE
  )
}

.bhf_flag_summary <- function(x, name) {
  values <- .bhf_draw_vector(x, name)
  if (any(!values %in% c(0, 1))) {
    .bhf_result_abort(
      "bhf_draw_contract_error",
      paste0("Posterior flag '", name, "' must contain only zero or one."),
      field = name
    )
  }
  logical_values <- as.logical(values)
  list(
    all = all(logical_values),
    any = any(logical_values),
    rate = mean(logical_values),
    draws = length(logical_values)
  )
}

.bhf_a_star_available <- function(fit) {
  identical(as.integer(fit$data$stan_data$use_deattenuation), 1L)
}

.bhf_sampling_variance_provenance <- function(fit) {
  if (is.list(fit$data$provenance) &&
      is.list(fit$data$provenance$sampling_variances)) {
    return(fit$data$provenance$sampling_variances)
  }
  if (is.list(fit$data$sampling_variance_info)) {
    return(fit$data$sampling_variance_info$provenance)
  }
  NULL
}

#' Extract Article-Aligned Variance Decomposition Results
#'
#' Summarizes separate latent diagnostics and probability-scale Estimands A,
#' A*, and B. No legacy 0.3.0 quantity is silently reinterpreted.
#'
#' @param fit A versioned `bhf_fit` object with schema `0.5.0`.
#' @param prob Credible-interval probability.
#' @param print Retained transition argument. If `TRUE`, prints the tidy table.
#' @return A `bhf_variance_decomposition` list with `latent`, `A`, `A_star`,
#'   `B`, `gaps`, and `summary_table` components. `A_star$available` is `FALSE`
#'   when de-attenuation was disabled.
#' @details Estimands A, A*, and B are probability-scale decompositions;
#'   latent-scale SDs, variances, and ICCs are reported separately. A* treats
#'   sampling variances as fixed inputs and does not propagate their estimation
#'   uncertainty. Intervals are pseudo-posterior credible intervals. Repeated-
#'   sample coverage has been measured only in a reduced balanced synthetic
#'   design and is not guaranteed more generally.
#' @export
variance_decomposition <- function(fit, prob = 0.95, print = FALSE) {
  interval <- .bhf_interval_contract(prob)
  if (!is.logical(print) || length(print) != 1L || is.na(print)) {
    .bhf_result_abort("bhf_argument_error",
                      "'print' must be one non-missing logical value.")
  }
  .bhf_validate_result_fit(fit)

  latent_map <- c(
    sd_state = "sigma_state",
    sd_stratum = "sigma_stratum",
    sd_psu = "sigma_psu",
    variance_state = "var_state_latent",
    variance_stratum = "var_stratum_latent",
    variance_psu = "var_psu_latent",
    variance_level1 = "var_level1_latent",
    variance_total = "var_total_latent",
    icc_state = "icc_state_latent",
    icc_stratum = "icc_stratum_latent",
    icc_psu = "icc_psu_latent"
  )
  a_map <- c(
    mean_probability = "p_bar_A", between = "var_between_A",
    within = "var_within_A", total = "var_total_A",
    proportion = "prop_between_A"
  )
  a_star_map <- c(
    mean_probability = "p_bar_A_star", between = "var_between_A_star",
    within = "var_within_A_star", total = "var_total_A_star",
    proportion = "prop_between_A_star"
  )
  b_map <- c(
    mean_probability = "p_bar_B", between = "var_between_B",
    within_binomial = "var_within_binomial_B",
    within_mixture = "var_within_mixture_B", within = "var_within_B",
    total = "var_total_B", proportion = "prop_between_B"
  )
  design_gap_map <- c(
    mean_probability = "gap_B_minus_A_mean",
    between = "gap_B_minus_A_between",
    within = "gap_B_minus_A_within",
    total = "gap_B_minus_A_total",
    proportion = "gap_B_minus_A_proportion"
  )
  noise_gap_map <- c(
    between = "gap_A_minus_A_star_between",
    total = "gap_A_minus_A_star_total",
    proportion = "gap_A_minus_A_star_proportion"
  )
  flags <- c(
    "A_proportion_defined", "A_star_at_boundary", "A_star_truncated",
    "A_star_proportion_defined", "B_proportion_defined",
    "mean_vhat_state"
  )
  required <- unique(c(
    unname(latent_map), unname(a_map), unname(a_star_map), unname(b_map),
    unname(design_gap_map), unname(noise_gap_map), flags
  ))
  draws <- .bhf_extract_required(fit, required)

  scalar_draws <- lapply(required, function(name) {
    .bhf_draw_vector(draws[[name]], name)
  })
  draw_counts <- unique(lengths(scalar_draws))
  if (length(draw_counts) != 1L) {
    .bhf_result_abort(
      "bhf_draw_shape_error",
      "Variance-decomposition posterior fields have inconsistent draw counts."
    )
  }

  latent_table <- .bhf_component_table(
    draws, latent_map, interval, "latent", "latent", "latent", TRUE
  )
  a_table <- .bhf_component_table(
    draws, a_map, interval, "estimand", "A", "probability", TRUE
  )
  b_table <- .bhf_component_table(
    draws, b_map, interval, "estimand", "B", "probability", TRUE
  )
  design_gap <- .bhf_component_table(
    draws, design_gap_map, interval, "gap", "B_minus_A", "probability", TRUE
  )

  a_star_available <- .bhf_a_star_available(fit)
  if (a_star_available) {
    a_star_table <- .bhf_component_table(
      draws, a_star_map, interval, "estimand", "A_star", "probability", TRUE
    )
    noise_gap <- .bhf_component_table(
      draws, noise_gap_map, interval, "gap", "A_minus_A_star",
      "probability", TRUE
    )
    correction <- as.list(.bhf_draw_summary(
      .bhf_draw_vector(draws$mean_vhat_state, "mean_vhat_state"),
      interval, "mean_vhat_state"
    ))
    a_star_flags <- list(
      at_boundary = .bhf_flag_summary(
        draws$A_star_at_boundary, "A_star_at_boundary"
      ),
      truncated = .bhf_flag_summary(
        draws$A_star_truncated, "A_star_truncated"
      ),
      proportion_defined = .bhf_flag_summary(
        draws$A_star_proportion_defined, "A_star_proportion_defined"
      )
    )
  } else {
    a_star_table <- NULL
    noise_gap <- NULL
    correction <- NULL
    a_star_flags <- NULL
  }

  a_star_summary_rows <- if (a_star_available) a_star_table else
    .bhf_unavailable_table(
      names(a_star_map), interval, "estimand", "A_star", "probability"
    )
  noise_gap_rows <- if (a_star_available) noise_gap else
    .bhf_unavailable_table(
      names(noise_gap_map), interval, "gap", "A_minus_A_star", "probability"
    )
  gap_table <- rbind(design_gap, noise_gap_rows)
  summary_table <- rbind(
    latent_table, a_table, a_star_summary_rows, b_table, gap_table
  )
  rownames(summary_table) <- NULL

  result <- list(
    schema_version = "0.5.0",
    interval = interval,
    latent = list(scale = "latent", summary = latent_table),
    A = list(
      available = TRUE, estimand = "A", scale = "probability",
      summary = a_table,
      proportion_defined = .bhf_flag_summary(
        draws$A_proportion_defined, "A_proportion_defined"
      )
    ),
    A_star = list(
      available = a_star_available,
      estimand = "A_star",
      scale = "probability",
      summary = a_star_table,
      correction = correction,
      flags = a_star_flags,
      provenance = .bhf_sampling_variance_provenance(fit),
      unavailable_reason = if (a_star_available) NULL else
        "deattenuation='none'; A* was not computed"
    ),
    B = list(
      available = TRUE, estimand = "B", scale = "probability",
      summary = b_table,
      proportion_defined = .bhf_flag_summary(
        draws$B_proportion_defined, "B_proportion_defined"
      )
    ),
    gaps = list(
      B_minus_A = design_gap,
      A_minus_A_star = noise_gap,
      table = gap_table
    ),
    summary_table = summary_table
  )
  class(result) <- c("bhf_variance_decomposition", "list")
  if (print) base::print(summary_table)
  invisible(result)
}

.bhf_domain_contract <- function(fit) {
  mapping <- fit$data$mapping$domain
  if (!is.data.frame(mapping) ||
      !all(c("label", "id") %in% names(mapping)) || !nrow(mapping) ||
      !identical(as.integer(mapping$id), seq_len(nrow(mapping)))) {
    .bhf_result_abort(
      "bhf_fit_schema_error",
      "The fit has an invalid domain label-ID mapping."
    )
  }
  shares <- as.numeric(fit$data$stan_data$w_state_pop_share)
  if (length(shares) != nrow(mapping) || anyNA(shares) ||
      any(!is.finite(shares)) || any(shares <= 0)) {
    .bhf_result_abort(
      "bhf_fit_schema_error",
      "The fit has invalid population shares for its domain mapping."
    )
  }
  n <- NULL
  domain_summary <- fit$data$domain_summary
  if (is.data.frame(domain_summary) &&
      all(c("state_id", "n") %in% names(domain_summary)) &&
      identical(as.integer(domain_summary$state_id), seq_len(nrow(mapping)))) {
    n <- as.integer(domain_summary$n)
  }
  if (is.null(n)) {
    state_id <- as.integer(fit$data$stan_data$state_id)
    n <- tabulate(state_id, nbins = nrow(mapping))
  }
  list(mapping = mapping, shares = shares, n = n)
}

.bhf_resolve_domain_estimand <- function(estimand, estimand_missing, type) {
  if (!is.null(type)) {
    if (!is.character(type) || length(type) != 1L || is.na(type) ||
        !type %in% c("conditional", "marginal")) {
      .bhf_result_abort(
        "bhf_deprecated_argument_error",
        "Deprecated 'type' must be exactly 'conditional' or 'marginal'.",
        argument = "type"
      )
    }
    if (!estimand_missing) {
      .bhf_result_abort(
        "bhf_argument_error",
        "Supply only 'estimand'; it cannot be combined with deprecated 'type'."
      )
    }
    if (identical(type, "marginal")) {
      .bhf_result_abort(
        "bhf_legacy_marginal_error",
        paste0(
          "Deprecated type='marginal' is not equivalent to Estimand B. ",
          "Refit with schema 0.5.0 and request estimand='B' explicitly."
        ),
        argument = "type"
      )
    }
    .bhf_result_warn(
      "bhf_deprecated_argument_warning",
      paste0(
        "type='conditional' is deprecated and maps to estimand='A'; ",
        "use estimand='A' before bhfvar 0.5.0."
      ),
      argument = "type"
    )
    return("A")
  }
  .bhf_match_choice(estimand, c("A", "B"), "estimand")
}

#' Extract Domain-Level Estimand A or B Probabilities
#'
#' @param fit A versioned `bhf_fit` object.
#' @param estimand Either `"A"` or `"B"`.
#' @param prob Credible-interval probability.
#' @param type Deprecated transition argument. `"conditional"` maps to A with
#'   a classed warning; `"marginal"` fails because it is not Estimand B.
#' @return A `bhf_domain_estimates` data frame in canonical domain-ID order,
#'   with label/ID, population share, sample size, posterior summaries, dynamic
#'   interval bounds, and estimand metadata.
#' @export
domain_estimates <- function(fit, estimand = c("A", "B"), prob = 0.95,
                             type = NULL) {
  estimand_missing <- missing(estimand)
  estimand <- .bhf_resolve_domain_estimand(estimand, estimand_missing, type)
  interval <- .bhf_interval_contract(prob)
  .bhf_validate_result_fit(fit)
  domain <- .bhf_domain_contract(fit)
  S <- nrow(domain$mapping)
  field <- if (identical(estimand, "A")) "p_state_A" else "p_state_B"
  draws <- .bhf_extract_required(fit, field)
  probabilities <- .bhf_draw_matrix(draws[[field]], S, field)

  summaries <- lapply(seq_len(S), function(s) {
    .bhf_draw_summary(probabilities[, s], interval, paste0(field, "[", s, "]"))
  })
  summary_matrix <- do.call(rbind, summaries)
  result <- data.frame(
    domain = as.character(domain$mapping$label),
    domain_id = as.integer(domain$mapping$id),
    estimand = rep(estimand, S),
    mean = summary_matrix[, "mean"],
    sd = summary_matrix[, "sd"],
    lower = summary_matrix[, "lower"],
    median = summary_matrix[, "median"],
    upper = summary_matrix[, "upper"],
    prob = rep(interval$prob, S),
    pop_share = domain$shares,
    n = domain$n,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  attr(result, "schema_version") <- "0.5.0"
  attr(result, "estimand") <- estimand
  attr(result, "interval") <- interval
  class(result) <- c("bhf_domain_estimates", "data.frame")
  result
}

#' Extract an Overall Estimand A or B Population Probability
#'
#' @param fit A versioned `bhf_fit` object.
#' @param estimand Either `"A"` or `"B"`.
#' @param prob Credible-interval probability.
#' @return A `bhf_overall_estimate` list containing the estimand, interval
#'   contract, posterior summary, and population-share provenance.
#' @export
overall_estimate <- function(fit, estimand = c("A", "B"), prob = 0.95) {
  estimand <- .bhf_match_choice(estimand, c("A", "B"), "estimand")
  interval <- .bhf_interval_contract(prob)
  .bhf_validate_result_fit(fit)
  domain <- .bhf_domain_contract(fit)
  suffix <- if (identical(estimand, "A")) "A" else "B"
  mean_field <- paste0("p_bar_", suffix)
  state_field <- paste0("p_state_", suffix)
  draws <- .bhf_extract_required(fit, c(mean_field, state_field))
  overall <- .bhf_draw_vector(draws[[mean_field]], mean_field)
  state <- .bhf_draw_matrix(draws[[state_field]], nrow(domain$mapping),
                            state_field)
  if (nrow(state) != length(overall)) {
    .bhf_result_abort(
      "bhf_draw_shape_error",
      "Overall and state-probability fields have inconsistent draw counts."
    )
  }
  reconstructed <- as.vector(state %*% domain$shares)
  if (!.bhf_close(overall, reconstructed)) {
    .bhf_result_abort(
      "bhf_result_invariant_error",
      paste0("Posterior field '", mean_field,
             "' does not equal the population-share reconstruction."),
      invariant = "overall.weighted_reconstruction"
    )
  }
  summary <- as.list(.bhf_draw_summary(overall, interval, mean_field))
  share_provenance <- fit$data$provenance$population_shares
  if (!is.list(share_provenance)) {
    .bhf_result_abort(
      "bhf_fit_schema_error",
      "The fit is missing population-share provenance."
    )
  }
  population_shares <- list(
    values = stats::setNames(
      domain$shares, as.character(domain$mapping$label)
    ),
    provenance = share_provenance
  )
  result <- c(
    list(
      schema_version = "0.5.0", estimand = estimand,
      scale = "probability", interval = interval,
      population_shares = population_shares
    ),
    summary
  )
  class(result) <- c("bhf_overall_estimate", "list")
  result
}

.bhf_log_lik_groups <- function(fit, aggregate, N) {
  analysis <- fit$data$analysis_data
  if (!is.data.frame(analysis) || nrow(analysis) != N) {
    .bhf_result_abort(
      "bhf_fit_schema_error",
      "Log-likelihood aggregation requires aligned analysis_data rows."
    )
  }
  if (identical(aggregate, "observation")) {
    labels <- if ("original_row" %in% names(analysis)) {
      as.character(analysis$original_row)
    } else as.character(seq_len(N))
    return(list(id = seq_len(N), labels = labels))
  }
  if (identical(aggregate, "psu")) {
    id <- as.integer(analysis$psu_flat_id)
    mapping <- fit$data$mapping$psu
    if (!is.data.frame(mapping) ||
        !all(c("id", "label", "stratum_label") %in% names(mapping))) {
      .bhf_result_abort("bhf_fit_schema_error",
                        "The fit has no valid PSU aggregation mapping.")
    }
    labels <- paste(mapping$stratum_label, mapping$label, sep = "::")
  } else {
    id <- as.integer(analysis$stratum_id)
    mapping <- fit$data$mapping$stratum
    if (!is.data.frame(mapping) ||
        !all(c("id", "label") %in% names(mapping))) {
      .bhf_result_abort("bhf_fit_schema_error",
                        "The fit has no valid stratum aggregation mapping.")
    }
    labels <- as.character(mapping$label)
  }
  if (length(id) != N || anyNA(id) || any(id < 1L) ||
      !identical(sort(unique(id)), seq_along(labels))) {
    .bhf_result_abort(
      "bhf_fit_schema_error",
      paste0("The ", aggregate, " aggregation IDs are not complete and aligned.")
    )
  }
  list(id = id, labels = labels)
}

.bhf_aggregate_log_lik <- function(log_lik, groups) {
  if (length(groups$id) != ncol(log_lik)) {
    .bhf_result_abort(
      "bhf_draw_shape_error",
      "Log-likelihood columns do not align with aggregation IDs."
    )
  }
  columns <- lapply(seq_along(groups$labels), function(group) {
    rowSums(log_lik[, groups$id == group, drop = FALSE])
  })
  result <- do.call(cbind, columns)
  colnames(result) <- groups$labels
  rownames(result) <- rownames(log_lik)
  result
}

#' Extract Raw or Pseudo Pointwise Log-Likelihood Contributions
#'
#' Aggregation only sums stored observation contributions. It does not refit
#' the model and does not establish ordinary observation- or cluster-level LOO.
#'
#' @param fit A versioned `bhf_fit` object.
#' @param kind `"pseudo"` (default) or `"raw"`.
#' @param aggregate `"observation"`, `"psu"`, or `"stratum"`.
#' @return A draw-by-unit `bhf_log_lik` matrix with explicit `kind`,
#'   `aggregate`, and interpretation metadata.
#' @details The default pseudo contributions target the fitted
#'   pseudo-likelihood. Aggregation is arithmetic only; this function does not
#'   establish ordinary observation- or cluster-level LOO validity.
#' @export
log_lik <- function(fit, kind = c("pseudo", "raw"),
                    aggregate = c("observation", "psu", "stratum")) {
  kind <- .bhf_match_choice(kind, c("pseudo", "raw"), "kind")
  aggregate <- .bhf_match_choice(
    aggregate, c("observation", "psu", "stratum"), "aggregate"
  )
  .bhf_validate_result_fit(fit)
  field <- if (identical(kind, "pseudo")) "log_lik_pseudo" else
    "log_lik_raw"
  draws <- .bhf_extract_required(fit, field)
  N <- as.integer(fit$data$stan_data$N)
  values <- .bhf_draw_matrix(draws[[field]], N, field)
  groups <- .bhf_log_lik_groups(fit, aggregate, N)
  result <- .bhf_aggregate_log_lik(values, groups)

  scope <- list(
    fitted_target = "survey-weighted Bayesian pseudo-likelihood",
    contribution = if (identical(kind, "pseudo"))
      "survey-weighted pseudo-log-likelihood" else
      "unweighted Bernoulli log-likelihood diagnostic",
    aggregation = aggregate,
    ordinary_loo_supported = FALSE,
    cluster_loo_supported = FALSE,
    caveat = paste0(
      "Aggregation sums contributions only; it is not leave-one-unit refitting ",
      "and does not provide standard LOO/WAIC guarantees for a pseudo-posterior."
    )
  )
  attr(result, "schema_version") <- "0.5.0"
  attr(result, "kind") <- kind
  attr(result, "aggregation") <- aggregate
  attr(result, "scope") <- scope
  class(result) <- c("bhf_log_lik", "matrix", "array")

  if (identical(kind, "pseudo")) {
    .bhf_result_warn(
      "bhf_pseudo_loo_warning",
      paste0(
        "Pseudo log-likelihood is survey-weighted and is not an ordinary ",
        "likelihood; standard LOO/WAIC interpretation is unsupported."
      ),
      scope = scope
    )
  } else {
    .bhf_result_warn(
      "bhf_raw_log_lik_warning",
      paste0(
        "Raw log-likelihood is unweighted and does not match the fitted ",
        "pseudo-likelihood target; standard LOO interpretation is unsupported."
      ),
      scope = scope
    )
  }
  result
}
