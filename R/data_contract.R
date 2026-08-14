#' Resolve Population Shares for the Observed Domain Universe
#'
#' Internal helper that validates externally supplied population shares or
#' estimates them from retained raw survey weights. Returned shares are named,
#' aligned to `domain_labels`, strictly positive, and normalized to sum to one.
#'
#' @param population_shares NULL or a named numeric vector.
#' @param domain_labels Observed domain labels in canonical state-ID order.
#' @param state_id Integer state IDs aligned with `raw_weights`.
#' @param raw_weights Positive finite raw survey weights.
#'
#' @return An internal `bhf_population_shares` object.
#' @keywords internal
resolve_population_shares <- function(population_shares,
                                      domain_labels,
                                      state_id,
                                      raw_weights) {
  if (!is.atomic(domain_labels) || length(domain_labels) == 0L ||
      anyNA(domain_labels)) {
    stop(
      "`domain_labels` must contain at least one non-missing observed domain.",
      call. = FALSE
    )
  }

  labels <- as.character(domain_labels)
  if (any(!nzchar(trimws(labels)))) {
    stop("`domain_labels` must contain non-blank labels.", call. = FALSE)
  }
  if (anyDuplicated(labels)) {
    stop("`domain_labels` must be unique.", call. = FALSE)
  }

  if (!is.numeric(state_id) || length(state_id) != length(raw_weights) ||
      anyNA(state_id) || any(!is.finite(state_id)) ||
      any(state_id != floor(state_id))) {
    stop(
      "`state_id` must be a finite integer-valued vector with one entry per raw weight.",
      call. = FALSE
    )
  }
  state_id <- as.integer(state_id)

  if (!is.numeric(raw_weights) || length(raw_weights) != length(state_id)) {
    stop(
      "`raw_weights` must be a numeric vector with one value per state ID.",
      call. = FALSE
    )
  }
  if (anyNA(raw_weights) || any(!is.finite(raw_weights)) ||
      any(raw_weights <= 0)) {
    stop(
      "`raw_weights` must contain only finite, strictly positive values.",
      call. = FALSE
    )
  }

  S <- length(labels)
  if (any(state_id < 1L | state_id > S)) {
    stop(
      "`state_id` values must be between 1 and `length(domain_labels)`.",
      call. = FALSE
    )
  }
  if (!identical(sort(unique(state_id)), seq_len(S))) {
    stop("`domain_labels` must contain only observed domains.", call. = FALSE)
  }

  if (is.null(population_shares)) {
    input_values <- vapply(
      seq_len(S),
      function(s) sum(raw_weights[state_id == s]),
      numeric(1)
    )
    source <- "estimated_from_raw_weights"
  } else {
    if (!is.numeric(population_shares) || !is.atomic(population_shares)) {
      stop("`population_shares` must be a numeric vector.", call. = FALSE)
    }

    share_names <- names(population_shares)
    if (is.null(share_names)) {
      stop("`population_shares` must be a named vector.", call. = FALSE)
    }
    if (anyNA(share_names) || any(!nzchar(trimws(share_names)))) {
      stop(
        "`population_shares` names must be non-missing and non-blank.",
        call. = FALSE
      )
    }
    if (anyDuplicated(share_names)) {
      stop("`population_shares` names must be unique.", call. = FALSE)
    }

    missing_names <- setdiff(labels, share_names)
    unknown_names <- setdiff(share_names, labels)
    if (length(missing_names) > 0L || length(unknown_names) > 0L) {
      details <- c(
        if (length(missing_names) > 0L) {
          paste0("Missing: ", paste(missing_names, collapse = ", "), ".")
        },
        if (length(unknown_names) > 0L) {
          paste0("Unknown: ", paste(unknown_names, collapse = ", "), ".")
        }
      )
      stop(
        "`population_shares` names must exactly match observed domain labels. ",
        paste(details, collapse = " "),
        call. = FALSE
      )
    }

    input_values <- as.numeric(population_shares[labels])
    if (anyNA(input_values) || any(!is.finite(input_values)) ||
        any(input_values <= 0)) {
      stop(
        "`population_shares` must contain only finite, strictly positive values.",
        call. = FALSE
      )
    }
    source <- "external_known"
  }

  input_sum <- sum(input_values)
  if (!is.finite(input_sum) || input_sum <= 0) {
    stop(
      "The population-share normalization sum must be finite and strictly positive.",
      call. = FALSE
    )
  }

  normalized_values <- input_values / input_sum
  if (anyNA(normalized_values) || any(!is.finite(normalized_values)) ||
      any(normalized_values <= 0)) {
    stop(
      "Normalized population shares must remain finite and strictly positive.",
      call. = FALSE
    )
  }

  values <- stats::setNames(normalized_values, labels)
  structure(
    list(
      values = values,
      source = source,
      input_sum = input_sum,
      normalization_factor = 1 / input_sum
    ),
    class = c("bhf_population_shares", "list")
  )
}

# Stable identifiers for the prepared-data contract. The package schema
# is a semantic version string, while the Stan-data schema is an integer that
# can later be required directly by the Stan program.
bhf_data_contract_spec <- function() {
  list(
    schema_version = "0.5.0",
    contract_id = "bhfvar-data-contract-0.5.0",
    stan_data_schema_version = 1L
  )
}

# Closed article prior-sensitivity catalog. Integer codes are the only values
# transported to Stan; the human-readable variant and hyperparameters remain
# in prepared-data provenance.
sigma_state_prior_catalog <- function() {
  list(
    half_t3_2.5 = list(
      variant = "half_t3_2.5", code = 1L, family = "half_student_t",
      df = 3, location = 0, scale = 2.5,
      article_label = "half-t(3, 0, 2.5)"
    ),
    half_normal_1 = list(
      variant = "half_normal_1", code = 2L, family = "half_normal",
      df = NULL, location = 0, scale = 1,
      article_label = "half-Normal(0, 1)"
    ),
    half_cauchy_2.5 = list(
      variant = "half_cauchy_2.5", code = 3L, family = "half_cauchy",
      df = NULL, location = 0, scale = 2.5,
      article_label = "half-Cauchy(0, 2.5)"
    ),
    half_t3_5 = list(
      variant = "half_t3_5", code = 4L, family = "half_student_t",
      df = 3, location = 0, scale = 5,
      article_label = "half-t(3, 0, 5)"
    )
  )
}

resolve_sigma_state_prior <- function(
    sigma_state_prior = names(sigma_state_prior_catalog())) {
  catalog <- sigma_state_prior_catalog()
  choices <- names(catalog)
  if (identical(sigma_state_prior, choices)) {
    sigma_state_prior <- choices[[1L]]
  }
  if (!is.character(sigma_state_prior) ||
      length(sigma_state_prior) != 1L || is.na(sigma_state_prior) ||
      !sigma_state_prior %in% choices) {
    stop(
      "sigma_state_prior must be exactly one of: ",
      paste(choices, collapse = ", "),
      call. = FALSE
    )
  }

  result <- catalog[[sigma_state_prior]]
  result$source <- if (identical(result$code, 1L)) {
    "default_baseline"
  } else {
    "article_sensitivity_selector"
  }
  result$baseline <- identical(result$code, 1L)
  result$varied_component <- "state"
  result$other_random_effect_priors_fixed <- TRUE
  result
}

# Assemble row provenance without enforcing cross-field invariants. Those
# checks belong to the dedicated data-invariant gate.
assemble_row_provenance <- function(n_original,
                                    n_used,
                                    retained_rows,
                                    dropped_rows,
                                    missing_reason_ledger) {
  list(
    n_original = as.integer(n_original),
    n_used = as.integer(n_used),
    retained_rows = as.integer(retained_rows),
    dropped_rows = as.integer(dropped_rows),
    missing_reason_ledger = missing_reason_ledger
  )
}

# Consolidate the independently produced provenance records into one stable
# surface while their existing transition fields remain available.
assemble_data_provenance <- function(contract,
                                     rows,
                                     weights,
                                     population_shares,
                                     sampling_variances,
                                     prior) {
  list(
    contract = contract,
    rows = rows,
    weights = weights,
    population_shares = population_shares,
    sampling_variances = sampling_variances,
    prior = prior
  )
}

.bhf_contract_abort <- function(invariant, message) {
  stop(structure(
    list(message = message, call = NULL, invariant = invariant),
    class = c("bhf_data_contract_error", "error", "condition")
  ))
}

.bhf_contract_require <- function(ok, invariant, message) {
  if (!isTRUE(ok)) {
    .bhf_contract_abort(invariant, message)
  }
  invisible(NULL)
}

.bhf_contract_close <- function(x, y) {
  length(x) == length(y) && all(
    is.finite(x) & is.finite(y) &
      abs(x - y) <= 1e-10 + 1e-8 * pmax(abs(x), abs(y))
  )
}

#' Validate the Complete Prepared-Data Contract
#'
#' Internal final gate for cross-field transport, mapping, and provenance
#' invariants. Field-local Stan support is delegated to `validate_stan_data()`.
#'
#' @param x A prepared `bhf_data` object.
#' @return `x`, invisibly, or a classed `bhf_data_contract_error`.
#' @keywords internal
validate_bhf_data_contract <- function(x) {
  assert_contract <- .bhf_contract_require
  spec <- bhf_data_contract_spec()
  required_top <- c(
    "schema_version", "contract_id", "stan_data_schema_version",
    "stan_data", "mapping", "row_provenance", "provenance",
    "prior_info", "weight_info", "population_share_info",
    "sampling_variance_info", "analysis_data", "domain_summary",
    "input_info"
  )

  assert_contract(inherits(x, "bhf_data"), "schema.class",
          "Prepared data must inherit from 'bhf_data'.")
  assert_contract(all(required_top %in% names(x)), "schema.fields",
          "Prepared data is missing required contract fields.")
  assert_contract(
    identical(x$schema_version, spec$schema_version) &&
      identical(x$contract_id, spec$contract_id) &&
      identical(x$stan_data_schema_version, spec$stan_data_schema_version) &&
      identical(x$stan_data$data_schema_version, spec$stan_data_schema_version) &&
      identical(x$input_info$schema_version, spec$schema_version) &&
      identical(x$input_info$contract_id, spec$contract_id) &&
      identical(x$input_info$stan_data_schema_version,
                spec$stan_data_schema_version) &&
      identical(x$provenance$contract, spec),
    "schema.version",
    "Prepared-data schema identifiers are inconsistent."
  )

  tryCatch(
    validate_stan_data(x$stan_data),
    error = function(error) {
      .bhf_contract_abort(
        "stan.local",
        paste0("Stan-data validation failed: ", conditionMessage(error))
      )
    }
  )

  stan <- x$stan_data
  analysis <- x$analysis_data
  N <- as.integer(stan$N)
  S <- as.integer(stan$S)
  H <- as.integer(stan$H)
  J <- as.integer(stan$J)
  analysis_fields <- c(
    "original_row", "y", "domain_label", "state_id", "stratum_label",
    "stratum_id", "psu_label", "psu_in_stratum_id", "psu_flat_id",
    "raw_weight", "w_lik"
  )
  assert_contract(is.data.frame(analysis) && all(analysis_fields %in% names(analysis)) &&
            nrow(analysis) == N && x$input_info$n_used == N &&
            x$row_provenance$n_used == N,
          "analysis.shape", "Analysis rows do not match Stan N.")
  assert_contract(
    identical(as.integer(analysis$y), as.integer(stan$y)) &&
      all(analysis$y %in% c(0L, 1L)),
    "outcome.transport", "Outcome values changed before Stan transport."
  )

  rows <- x$row_provenance
  retained <- rows$retained_rows
  dropped <- rows$dropped_rows
  assert_contract(
    identical(x$provenance$rows, rows) &&
      identical(x$input_info$retained_rows, retained) &&
      identical(x$input_info$dropped_rows, dropped) &&
      identical(x$weight_info$original_row, retained) &&
      identical(analysis$original_row, retained) &&
      identical(rows$n_original, as.integer(x$input_info$n_original)) &&
      identical(rows$n_used, as.integer(x$input_info$n_used)) &&
      identical(retained, sort(unique(retained))) &&
      identical(dropped, sort(unique(dropped))) &&
      length(intersect(retained, dropped)) == 0L &&
      identical(sort(c(retained, dropped)), seq_len(rows$n_original)),
    "rows.partition", "Retained and dropped row provenance is inconsistent."
  )
  ledger <- rows$missing_reason_ledger
  role_columns <- c(
    outcome = x$input_info$outcome, domain = x$input_info$domain,
    strata = x$input_info$strata, psu = x$input_info$psu,
    weights = x$input_info$weights
  )
  ledger_ok <- is.data.frame(ledger) &&
    all(c("original_row", "field", "variable") %in% names(ledger)) &&
    all(ledger$original_row %in% dropped) &&
    all(dropped %in% ledger$original_row) &&
    all(ledger$field %in% names(role_columns)) &&
    all(ledger$variable == unname(role_columns[ledger$field])) &&
    !anyDuplicated(ledger[c("original_row", "field")])
  assert_contract(ledger_ok, "rows.ledger", "Missing-row ledger is inconsistent.")

  domain_map <- x$mapping$domain
  assert_contract(
    is.data.frame(domain_map) && nrow(domain_map) == S &&
      identical(domain_map$id, seq_len(S)) &&
      !anyDuplicated(domain_map$label) &&
      identical(as.integer(analysis$state_id), as.integer(stan$state_id)) &&
      identical(analysis$domain_label,
                as.character(domain_map$label[analysis$state_id])) &&
      identical(sort(unique(as.integer(analysis$state_id))), seq_len(S)),
    "mapping.domain", "Domain IDs and labels do not round-trip."
  )
  stratum_map <- x$mapping$stratum
  assert_contract(
    is.data.frame(stratum_map) && nrow(stratum_map) == H &&
      identical(stratum_map$id, seq_len(H)) &&
      !anyDuplicated(stratum_map$label) &&
      identical(as.integer(analysis$stratum_id), as.integer(stan$stratum_id)) &&
      identical(analysis$stratum_label,
                as.character(stratum_map$label[analysis$stratum_id])) &&
      identical(sort(unique(as.integer(analysis$stratum_id))), seq_len(H)),
    "mapping.stratum", "Stratum IDs and labels do not round-trip."
  )
  psu_map <- x$mapping$psu
  psu_fields <- c(
    "label", "id", "stratum_label", "stratum_id", "within_stratum_id"
  )
  mapped_psu <- if (is.data.frame(psu_map) && nrow(psu_map) == J &&
                    all(psu_fields %in% names(psu_map))) {
    psu_map[analysis$psu_flat_id, , drop = FALSE]
  } else {
    NULL
  }
  psu_ok <- !is.null(mapped_psu) &&
    identical(psu_map$id, seq_len(J)) &&
    !anyDuplicated(psu_map[c("stratum_id", "label")]) &&
    identical(psu_map$stratum_label,
              as.character(stratum_map$label[psu_map$stratum_id])) &&
    identical(as.integer(tabulate(psu_map$stratum_id, H)),
              as.integer(stan$J_h)) &&
    identical(analysis$psu_label, as.character(mapped_psu$label)) &&
    identical(as.integer(analysis$stratum_id),
              as.integer(mapped_psu$stratum_id)) &&
    identical(as.integer(analysis$psu_in_stratum_id),
              as.integer(mapped_psu$within_stratum_id)) &&
    identical(as.integer(analysis$psu_flat_id),
              as.integer(stan$psu_flat_id)) &&
    identical(sort(unique(as.integer(analysis$psu_flat_id))), seq_len(J))
  assert_contract(psu_ok, "mapping.psu", "Composite PSU IDs do not round-trip.")

  weights <- x$weight_info
  raw <- analysis$raw_weight
  likelihood <- analysis$w_lik
  assert_contract(
    identical(raw, weights$raw) &&
      identical(likelihood, weights$likelihood) &&
      identical(likelihood, as.vector(stan$w_lik)),
    "weights.transport", "Raw or likelihood weights changed between surfaces."
  )
  wp <- weights$provenance
  assert_contract(
    identical(x$provenance$weights, wp) &&
      identical(x$input_info$weight_scaling, wp$method) &&
      identical(wp$normalization, "global_sum_N_in_R") &&
      identical(as.integer(wp$N), N) &&
      .bhf_contract_close(wp$raw_sum, sum(raw)) &&
      .bhf_contract_close(wp$likelihood_sum, sum(likelihood)) &&
      abs(sum(likelihood) - N) <= 1e-10 * max(1, N),
    "weights.provenance", "Weight provenance is inconsistent."
  )
  if (identical(wp$method, "mean_one")) {
    expected <- (raw / max(raw)) / mean(raw / max(raw))
    assert_contract(is.null(wp$legacy_domain_factors) &&
              .bhf_contract_close(likelihood, expected) &&
              .bhf_contract_close(wp$global_factor, N / sum(raw)),
            "weights.normalization", "Mean-one weight normalization is invalid.")
  } else if (identical(wp$method, "legacy_d2")) {
    factors <- wp$legacy_domain_factors
    labels <- as.character(domain_map$label)
    expected_input <- raw * unname(factors[labels[analysis$state_id]])
    assert_contract(identical(names(factors), labels) && all(is.finite(factors)) &&
              all(factors > 0) &&
              .bhf_contract_close(likelihood, expected_input / mean(expected_input)) &&
              .bhf_contract_close(wp$global_factor, N / sum(expected_input)),
            "weights.normalization", "Legacy D2 weight normalization is invalid.")
  } else {
    .bhf_contract_abort("weights.provenance", "Unknown weight-scaling method.")
  }

  share_info <- x$population_share_info
  share_values <- share_info$values
  labels <- as.character(domain_map$label)
  assert_contract(inherits(share_info, "bhf_population_shares") &&
            identical(names(share_values), labels) &&
            all(is.finite(share_values)) && all(share_values > 0) &&
            abs(sum(share_values) - 1) <= 1e-12,
          "shares.labels", "Population-share labels or support are invalid.")
  assert_contract(.bhf_contract_close(unname(share_values),
                              as.vector(stan$w_state_pop_share)) &&
            .bhf_contract_close(unname(share_values), x$domain_summary$pop_share),
          "shares.transport", "Population shares changed between surfaces.")
  share_prov <- x$provenance$population_shares
  assert_contract(identical(share_prov$source, share_info$source) &&
            identical(share_prov$domain_labels, labels) &&
            identical(x$input_info$population_share_source, share_info$source) &&
            .bhf_contract_close(share_prov$input_sum, share_info$input_sum) &&
            .bhf_contract_close(share_info$normalization_factor,
                                1 / share_info$input_sum),
          "shares.provenance", "Population-share provenance is inconsistent.")
  if (identical(share_info$source, "estimated_from_raw_weights")) {
    totals <- vapply(seq_len(S), function(s) sum(raw[analysis$state_id == s]),
                     numeric(1))
    assert_contract(.bhf_contract_close(unname(share_values), totals / sum(totals)),
            "shares.provenance", "Estimated shares do not use raw weights.")
  } else {
    assert_contract(identical(share_info$source, "external_known"),
            "shares.provenance", "Unknown population-share source.")
  }

  vhat <- x$sampling_variance_info
  vp <- vhat$provenance
  assert_contract(identical(x$provenance$sampling_variances, vp) &&
            identical(vp$domain_labels, labels) &&
            identical(vp$units, "probability_squared") &&
            identical(vp$uncertainty_propagated, FALSE),
          "vhat.provenance", "Sampling-variance provenance is inconsistent.")
  assert_contract(identical(x$input_info$deattenuation, vp$mode) &&
            identical(x$input_info$use_deattenuation, vhat$enabled) &&
            identical(as.integer(stan$use_deattenuation),
                      as.integer(vhat$enabled)),
          "vhat.mode", "De-attenuation mode flags are inconsistent.")
  assert_contract(.bhf_contract_close(as.vector(stan$vhat_state), vhat$stan_values),
          "vhat.transport", "Sampling variances changed before Stan transport.")
  if (identical(vp$mode, "none")) {
    assert_contract(identical(vhat$enabled, FALSE) &&
              identical(vhat$a_star_available, FALSE) &&
              is.null(vhat$named_values) && all(vhat$stan_values == 0) &&
              all(is.na(x$domain_summary$vhat)) &&
              identical(vp$source, "none") &&
              identical(vp$fixed_input, FALSE) &&
              identical(vp$stan_placeholder, TRUE),
            "vhat.mode", "None-mode placeholder contract is invalid.")
  } else {
    assert_contract(vp$mode %in% c("supplied", "taylor") &&
              identical(vhat$enabled, TRUE) &&
              identical(vhat$a_star_available, TRUE) &&
              identical(names(vhat$named_values), labels) &&
              .bhf_contract_close(unname(vhat$named_values), vhat$stan_values) &&
              .bhf_contract_close(x$domain_summary$vhat, vhat$stan_values) &&
              identical(vp$fixed_input, TRUE) &&
              identical(vp$stan_placeholder, FALSE),
            "vhat.mode", "Enabled sampling-variance contract is invalid.")
    expected_source <- if (identical(vp$mode, "supplied")) {
      "external"
    } else {
      "builtin_survey_taylor"
    }
    assert_contract(identical(vp$source, expected_source),
            "vhat.provenance", "Sampling-variance source is invalid.")
  }

  prior <- x$prior_info
  assert_contract(identical(x$provenance$prior, prior) &&
            prior$alpha$mean_source %in%
              c("estimated_from_raw_weighted_outcome_logit", "user_supplied") &&
            prior$alpha$sd_source %in% c("default", "user_supplied") &&
            identical(prior$alpha$scale, "logit"),
          "prior.provenance", "Prior provenance is inconsistent.")
  assert_contract(.bhf_contract_close(prior$alpha$mean, stan$prior_alpha_mean) &&
            .bhf_contract_close(prior$alpha$sd, stan$prior_alpha_sd),
          "prior.transport", "Alpha prior changed before Stan transport.")
  state_prior <- prior$sigma_state
  expected_state_prior <- tryCatch(
    resolve_sigma_state_prior(state_prior$variant),
    error = function(error) NULL
  )
  assert_contract(
    !is.null(expected_state_prior) &&
      identical(state_prior, expected_state_prior) &&
      identical(x$input_info$sigma_state_prior, state_prior$variant) &&
      identical(as.integer(stan$sigma_state_prior_code), state_prior$code),
    "prior.state_selector",
    "The sigma_state prior selector or provenance is inconsistent."
  )
  prevalence <- sum(raw * analysis$y) / sum(raw)
  assert_contract(.bhf_contract_close(prior$alpha$design_weighted_prevalence,
                              prevalence) &&
            (prior$alpha$mean_source !=
               "estimated_from_raw_weighted_outcome_logit" ||
               .bhf_contract_close(prior$alpha$mean, stats::qlogis(prevalence))) &&
            (prior$alpha$sd_source != "default" ||
               .bhf_contract_close(prior$alpha$sd, 0.5)) &&
            identical(prior$random_effect_sd$family, "half_student_t") &&
            identical(prior$random_effect_sd$df, 3) &&
            identical(prior$random_effect_sd$location, 0) &&
            identical(prior$random_effect_sd$scale, 2.5) &&
            identical(prior$random_effect_sd$applies_to,
                      c("state", "stratum", "psu")),
          "prior.provenance", "Prior values do not match their provenance.")

  summary <- x$domain_summary
  counts <- tabulate(analysis$state_id, S)
  weight_sums <- vapply(seq_len(S), function(s) sum(raw[analysis$state_id == s]),
                        numeric(1))
  prop <- vapply(seq_len(S), function(s) {
    idx <- analysis$state_id == s
    sum(raw[idx] * analysis$y[idx]) / sum(raw[idx])
  }, numeric(1))
  eff_n <- vapply(seq_len(S), function(s) {
    w <- raw[analysis$state_id == s]
    sum(w)^2 / sum(w^2)
  }, numeric(1))
  assert_contract(is.data.frame(summary) && nrow(summary) == S &&
            identical(summary$state_id, seq_len(S)) &&
            identical(summary$n, counts) && identical(summary$n_obs, counts) &&
            sum(summary$n) == N &&
            .bhf_contract_close(summary$weight_sum, weight_sums) &&
            .bhf_contract_close(summary$prop_weighted, prop) &&
            .bhf_contract_close(summary$eff_n, eff_n),
          "summary.transport", "Domain summary does not match analysis rows.")

  invisible(x)
}
