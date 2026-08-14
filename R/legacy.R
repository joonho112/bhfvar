# Legacy-object detection and fail-closed migration policy.
#
# Version 0.3.0 used a two-random-effect model and assigned different meanings
# to A, A*, and B. Its prepared data and posterior draws therefore cannot be
# upgraded by renaming fields. Detection is supported for one release so that
# users receive a deterministic re-prepare/refit instruction.

.bhf_legacy_policy <- function() {
  list(
    legacy_schema_version = "0.3.0",
    current_schema_version = bhf_data_contract_spec()$schema_version,
    detection_support_release = "0.5.0",
    removal_eligible_release = "0.5.0",
    support_window = "one_minor_release",
    behavior = "detect_and_fail_closed",
    automatic_adapter = FALSE,
    silent_reinterpretation = FALSE,
    data_action = "reprepare",
    fit_action = "refit"
  )
}

.bhf_legacy_old_name_inventory <- function() {
  list(
    data_top_level = c(
      "stan_data", "mapping", "domain_summary", "input_info"
    ),
    stan_data = c(
      "N", "S", "H", "J", "y", "state_id", "stratum_id", "J_h",
      "psu_start", "psu_in_stratum_id", "w_lik", "w_state_pop_share",
      "prior_alpha_mean", "prior_alpha_sd", "use_deattenuation",
      "vhat_state"
    ),
    fit_top_level = c(
      "stanfit", "data", "model", "call", "diagnostics"
    ),
    posterior_draws = c(
      "alpha", "z_state", "z_psu", "sigma_state", "sigma_psu",
      "u_state", "u_psu", "eta", "var_between_state", "var_psu",
      "var_logistic", "var_within_state", "var_total_logit", "icc_state",
      "c_factor", "p_state_conditional", "p_state_marginal", "p_overall",
      "var_between_prob", "var_within_prob", "var_total_prob", "icc_prob",
      "var_between_deatten", "icc_deatten", "reliability_state",
      "reliability_avg", "log_lik", "y_rep"
    )
  )
}

.bhf_legacy_detection <- function(object_type,
                                  schema_version,
                                  status,
                                  contract_id = NA_character_,
                                  reason,
                                  inferred = FALSE) {
  structure(
    list(
      object_type = object_type,
      schema_version = schema_version,
      contract_id = contract_id,
      status = status,
      inferred = isTRUE(inferred),
      reason = reason
    ),
    class = c("bhf_schema_detection", "list")
  )
}

.bhf_legacy_scalar_character <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.bhf_legacy_has_exact_names <- function(x, required) {
  is.list(x) && all(required %in% names(x))
}

.bhf_legacy_is_data_signature <- function(x) {
  inventory <- .bhf_legacy_old_name_inventory()
  current_only <- c(
    "schema_version", "contract_id", "stan_data_schema_version",
    "analysis_data", "row_provenance", "provenance", "weight_info",
    "population_share_info", "sampling_variance_info", "prior_info"
  )
  stan_current_only <- c(
    "data_schema_version", "psu_flat_id", "sigma_state_prior_code"
  )

  inherits(x, "bhf_data") &&
    .bhf_legacy_has_exact_names(x, inventory$data_top_level) &&
    .bhf_legacy_has_exact_names(x$stan_data, inventory$stan_data) &&
    !any(current_only %in% names(x)) &&
    !any(stan_current_only %in% names(x$stan_data))
}

#' Detect the Schema of a BHF Object
#'
#' Classifies current, legacy 0.3.0, unsupported, and non-BHF objects without
#' mutating or adapting them. Legacy classification is inferred from the
#' frozen 0.3.0 structural signature because those objects predate explicit
#' schema markers.
#'
#' @param x An object to inspect.
#'
#' @return A `bhf_schema_detection` list with object type, schema version,
#'   status, contract identifier, inference flag, and reason.
#'
#' @keywords internal
detect_bhf_object_schema <- function(x) {
  spec <- bhf_data_contract_spec()
  fit_spec <- if (exists(".bhf_fit_contract_spec", mode = "function")) {
    .bhf_fit_contract_spec()
  } else {
    list(
      schema_version = spec$schema_version,
      contract_id = "bhfvar-fit-contract-0.5.0"
    )
  }

  if (inherits(x, "bhf_data")) {
    current_markers <-
      identical(x$schema_version, spec$schema_version) &&
      identical(x$contract_id, spec$contract_id) &&
      identical(x$stan_data_schema_version, spec$stan_data_schema_version) &&
      is.list(x$stan_data) &&
      identical(x$stan_data$data_schema_version,
                spec$stan_data_schema_version)

    if (current_markers) {
      return(.bhf_legacy_detection(
        "bhf_data", spec$schema_version, "current", spec$contract_id,
        "explicit current prepared-data markers"
      ))
    }
    if (.bhf_legacy_is_data_signature(x)) {
      return(.bhf_legacy_detection(
        "bhf_data", "0.3.0", "legacy", NA_character_,
        "inferred from the frozen pre-schema 0.3.0 field signature",
        inferred = TRUE
      ))
    }

    declared <- if (.bhf_legacy_scalar_character(x$schema_version)) {
      x$schema_version
    } else {
      NA_character_
    }
    contract <- if (.bhf_legacy_scalar_character(x$contract_id)) {
      x$contract_id
    } else {
      NA_character_
    }
    return(.bhf_legacy_detection(
      "bhf_data", declared, "unsupported", contract,
      "bhf_data markers do not match a supported current or legacy schema"
    ))
  }

  if (inherits(x, "bhf_fit")) {
    if (!is.list(x$data)) {
      return(.bhf_legacy_detection(
        "bhf_fit", NA_character_, "unsupported", NA_character_,
        "bhf_fit does not contain a detectable prepared-data object"
      ))
    }

    data_detection <- detect_bhf_object_schema(x$data)
    declared <- if (.bhf_legacy_scalar_character(x$schema_version)) {
      x$schema_version
    } else {
      data_detection$schema_version
    }
    legacy_marker_ok <- is.null(x$schema_version) && is.null(x$contract_id)
    current_marker_ok <-
      identical(x$schema_version, fit_spec$schema_version) &&
      identical(x$contract_id, fit_spec$contract_id)

    if (identical(data_detection$status, "legacy") && legacy_marker_ok) {
      return(.bhf_legacy_detection(
        "bhf_fit", "0.3.0", "legacy", NA_character_,
        "inferred from the nested frozen 0.3.0 bhf_data signature",
        inferred = TRUE
      ))
    }
    if (identical(data_detection$status, "current") && current_marker_ok) {
      return(.bhf_legacy_detection(
        "bhf_fit", fit_spec$schema_version, "current", fit_spec$contract_id,
        "fit and nested prepared data resolve to the current schema"
      ))
    }
    return(.bhf_legacy_detection(
      "bhf_fit", declared, "unsupported", NA_character_,
      "fit markers and nested prepared-data schema are inconsistent"
    ))
  }

  .bhf_legacy_detection(
    "unknown", NA_character_, "not_bhf", NA_character_,
    "object does not inherit from bhf_data or bhf_fit"
  )
}

.bhf_legacy_abort <- function(detection, expected_type) {
  policy <- .bhf_legacy_policy()
  legacy <- identical(detection$status, "legacy")
  action <- if (identical(expected_type, "bhf_data")) {
    "re-prepare from the original row-level data with prepare_bhf_data()"
  } else {
    paste0(
      "re-prepare from the original row-level data with prepare_bhf_data() ",
      "and refit the model with bhf_fit()"
    )
  }
  version <- if (is.na(detection$schema_version)) {
    "unknown"
  } else {
    detection$schema_version
  }
  message <- if (legacy) {
    paste0(
      "Legacy ", expected_type, " schema ", version,
      " cannot be used as schema ", policy$current_schema_version,
      " because the model and estimand semantics changed; ", action,
      ". Automatic adaptation and silent reinterpretation are unavailable."
    )
  } else {
    paste0(
      "Unsupported ", expected_type, " schema ", version, "; ", action,
      ". Automatic adaptation and silent reinterpretation are unavailable."
    )
  }

  condition <- structure(
    list(
      message = message,
      call = NULL,
      object_type = expected_type,
      detected_schema = detection$schema_version,
      current_schema = policy$current_schema_version,
      detected_status = detection$status,
      required_action = if (identical(expected_type, "bhf_data")) {
        policy$data_action
      } else {
        policy$fit_action
      },
      automatic_adapter = FALSE,
      support_release = policy$detection_support_release,
      removal_eligible_release = policy$removal_eligible_release
    ),
    class = c(
      if (legacy) "bhf_legacy_object_error" else
        "bhf_unsupported_schema_error",
      "bhf_schema_error", "error", "condition"
    )
  )
  stop(condition)
}

.bhf_legacy_assert_current_bhf_data <- function(x) {
  detection <- detect_bhf_object_schema(x)
  if (!identical(detection$object_type, "bhf_data") ||
      !identical(detection$status, "current")) {
    .bhf_legacy_abort(detection, "bhf_data")
  }
  validate_bhf_data_contract(x)
  invisible(x)
}

.bhf_legacy_assert_current_bhf_fit <- function(x) {
  detection <- detect_bhf_object_schema(x)
  if (!identical(detection$object_type, "bhf_fit") ||
      !identical(detection$status, "current")) {
    .bhf_legacy_abort(detection, "bhf_fit")
  }
  .bhf_legacy_assert_current_bhf_data(x$data)
  invisible(x)
}
