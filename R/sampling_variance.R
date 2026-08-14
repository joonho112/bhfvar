# Sampling-variance input contracts -----------------------------------------

.bhf_vhat_abort <- function(message, class) {
  condition <- structure(
    list(message = message, call = NULL),
    class = c(class, "error", "condition")
  )
  stop(condition)
}

.bhf_vhat_warn <- function(message, class) {
  condition <- structure(
    list(message = message, call = NULL),
    class = c(class, "warning", "condition")
  )
  warning(condition)
}

.bhf_survey_available <- function() {
  requireNamespace("survey", quietly = TRUE)
}

.validate_domain_labels <- function(domain_labels) {
  if (is.null(domain_labels) || length(domain_labels) == 0L ||
      !is.atomic(domain_labels)) {
    .bhf_vhat_abort(
      "domain_labels must be a non-empty atomic vector.",
      "bhf_vhat_alignment_error"
    )
  }

  labels <- as.character(domain_labels)
  if (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
    .bhf_vhat_abort(
      "domain_labels must be non-missing, non-blank, and unique.",
      "bhf_vhat_alignment_error"
    )
  }

  labels
}

.validate_supplied_vhat <- function(sampling_variances, domain_labels) {
  if (is.null(sampling_variances)) {
    .bhf_vhat_abort(
      "sampling_variances is required when deattenuation = 'supplied'.",
      "bhf_vhat_mode_error"
    )
  }
  if (!is.numeric(sampling_variances) || is.factor(sampling_variances) ||
      length(sampling_variances) == 0L) {
    .bhf_vhat_abort(
      "sampling_variances must be a non-empty numeric vector.",
      "bhf_vhat_value_error"
    )
  }
  if (any(!is.finite(sampling_variances)) || any(sampling_variances < 0)) {
    .bhf_vhat_abort(
      "sampling_variances must contain only finite, nonnegative values.",
      "bhf_vhat_value_error"
    )
  }

  supplied_names <- names(sampling_variances)
  if (is.null(supplied_names) || length(supplied_names) != length(sampling_variances) ||
      anyNA(supplied_names) || any(!nzchar(supplied_names)) ||
      anyDuplicated(supplied_names)) {
    .bhf_vhat_abort(
      "sampling_variances must have non-missing, non-blank, unique domain names.",
      "bhf_vhat_alignment_error"
    )
  }

  missing_labels <- setdiff(domain_labels, supplied_names)
  extra_labels <- setdiff(supplied_names, domain_labels)
  if (length(missing_labels) > 0L || length(extra_labels) > 0L) {
    details <- c(
      if (length(missing_labels) > 0L) {
        paste0("missing: ", paste(missing_labels, collapse = ", "))
      },
      if (length(extra_labels) > 0L) {
        paste0("unknown: ", paste(extra_labels, collapse = ", "))
      }
    )
    .bhf_vhat_abort(
      paste0(
        "sampling_variances names must exactly match the observed domain universe (",
        paste(details, collapse = "; "), ")."
      ),
      "bhf_vhat_alignment_error"
    )
  }

  aligned <- as.double(sampling_variances[match(domain_labels, supplied_names)])
  stats::setNames(aligned, domain_labels)
}

.validate_supplied_vhat_method <- function(sampling_variance_method) {
  supported <- c(
    "external_taylor",
    "external_replicate",
    "external_other"
  )

  if (is.null(sampling_variance_method) ||
      length(sampling_variance_method) != 1L ||
      !is.character(sampling_variance_method) ||
      is.na(sampling_variance_method) ||
      !sampling_variance_method %in% supported) {
    .bhf_vhat_abort(
      paste0(
        "sampling_variance_method must be one of: ",
        paste(supported, collapse = ", "),
        " when deattenuation = 'supplied'."
      ),
      "bhf_vhat_mode_error"
    )
  }

  sampling_variance_method
}

.validate_taylor_analysis_data <- function(analysis_data, domain_labels) {
  required <- c(
    "y",
    "domain_label",
    "stratum_id",
    "psu_flat_id",
    "raw_weight"
  )

  if (!is.data.frame(analysis_data) || nrow(analysis_data) == 0L) {
    .bhf_vhat_abort(
      "analysis_data must be a non-empty data frame for Taylor estimation.",
      "bhf_taylor_design_error"
    )
  }
  missing_columns <- setdiff(required, names(analysis_data))
  if (length(missing_columns) > 0L) {
    .bhf_vhat_abort(
      paste0(
        "analysis_data is missing required Taylor columns: ",
        paste(missing_columns, collapse = ", "),
        "."
      ),
      "bhf_taylor_design_error"
    )
  }

  result <- analysis_data[, required, drop = FALSE]
  if (!is.numeric(result$y) || anyNA(result$y) ||
      any(!is.finite(result$y)) || any(!result$y %in% c(0, 1))) {
    .bhf_vhat_abort(
      "analysis_data$y must contain only complete binary 0/1 values.",
      "bhf_taylor_design_error"
    )
  }

  observed_labels <- as.character(result$domain_label)
  if (anyNA(observed_labels) || any(!nzchar(observed_labels))) {
    .bhf_vhat_abort(
      "analysis_data$domain_label must be complete and non-blank.",
      "bhf_taylor_design_error"
    )
  }
  if (!setequal(unique(observed_labels), domain_labels)) {
    .bhf_vhat_abort(
      "analysis_data domains must exactly match domain_labels.",
      c("bhf_taylor_design_error", "bhf_vhat_alignment_error")
    )
  }

  integer_columns <- c("stratum_id", "psu_flat_id")
  for (column in integer_columns) {
    values <- result[[column]]
    if (!is.numeric(values) || anyNA(values) || any(!is.finite(values)) ||
        any(values < 1) || any(values != floor(values))) {
      .bhf_vhat_abort(
        paste0("analysis_data$", column, " must contain positive integer IDs."),
        "bhf_taylor_design_error"
      )
    }
    result[[column]] <- as.integer(values)
  }

  if (!identical(
    sort(unique(result$stratum_id)),
    seq_len(max(result$stratum_id))
  )) {
    .bhf_vhat_abort(
      "analysis_data$stratum_id must be consecutive from 1.",
      "bhf_taylor_design_error"
    )
  }
  psu_strata <- split(result$stratum_id, result$psu_flat_id)
  if (any(vapply(psu_strata, function(x) length(unique(x)) != 1L, logical(1)))) {
    .bhf_vhat_abort(
      "Each composite psu_flat_id must belong to exactly one stratum.",
      "bhf_taylor_design_error"
    )
  }

  psus_per_stratum <- vapply(
    split(result$psu_flat_id, result$stratum_id),
    function(x) length(unique(x)),
    integer(1)
  )
  if (any(psus_per_stratum < 2L)) {
    lonely <- names(psus_per_stratum)[psus_per_stratum < 2L]
    .bhf_vhat_abort(
      paste0(
        "Taylor estimation requires at least two PSUs per stratum; singleton ",
        "strata: ", paste(lonely, collapse = ", "), "."
      ),
      c("bhf_taylor_singleton_error", "bhf_taylor_design_error")
    )
  }

  if (!is.numeric(result$raw_weight) || anyNA(result$raw_weight) ||
      any(!is.finite(result$raw_weight)) || any(result$raw_weight <= 0)) {
    .bhf_vhat_abort(
      "analysis_data$raw_weight must contain finite, strictly positive values.",
      "bhf_taylor_design_error"
    )
  }

  result$domain_label <- factor(observed_labels, levels = domain_labels)
  result
}

#' Estimate Domain Sampling Variances by Taylor Linearization
#'
#' Uses a one-stage, stratified PSU survey design and raw survey weights to
#' estimate the design-based variance of each domain Hajek proportion. The
#' function deliberately fails for singleton strata and does not apply caps,
#' floors, or heuristic fallbacks.
#'
#' @keywords internal
estimate_taylor_vhat <- function(analysis_data, domain_labels) {
  if (!.bhf_survey_available()) {
    .bhf_vhat_abort(
      paste0(
        "Taylor sampling-variance estimation requires the suggested package ",
        "'survey'. Install it or use deattenuation = 'supplied'/'none'."
      ),
      "bhf_taylor_dependency_error"
    )
  }

  survey_option_names <- c(
    "survey.lonely.psu",
    "survey.adjust.domain.lonely"
  )
  all_previous_options <- options()
  option_was_set <- survey_option_names %in% names(all_previous_options)
  previous_values <- lapply(
    survey_option_names,
    function(name) all_previous_options[[name]]
  )
  restore_survey_options <- function() {
    for (index in seq_along(survey_option_names)) {
      value <- if (option_was_set[index]) previous_values[[index]] else NULL
      options(stats::setNames(list(value), survey_option_names[index]))
    }
  }
  on.exit(restore_survey_options(), add = TRUE)
  options(
    survey.lonely.psu = "fail",
    survey.adjust.domain.lonely = FALSE
  )

  labels <- .validate_domain_labels(domain_labels)
  survey_data <- .validate_taylor_analysis_data(analysis_data, labels)

  estimates <- tryCatch(
    {
      design <- survey::svydesign(
        ids = ~psu_flat_id,
        strata = ~stratum_id,
        weights = ~raw_weight,
        data = survey_data,
        nest = TRUE
      )
      survey::svyby(
        ~y,
        ~domain_label,
        design,
        survey::svymean,
        vartype = "se",
        keep.var = TRUE,
        covmat = TRUE,
        na.rm = FALSE,
        drop.empty.groups = FALSE
      )
    },
    error = function(error) {
      .bhf_vhat_abort(
        paste0("Taylor sampling-variance estimation failed: ", conditionMessage(error)),
        "bhf_taylor_design_error"
      )
    }
  )

  result_labels <- as.character(estimates$domain_label)
  if (anyNA(result_labels) || anyDuplicated(result_labels) ||
      !setequal(result_labels, labels)) {
    .bhf_vhat_abort(
      "Taylor output labels do not match the observed domain universe.",
      "bhf_vhat_alignment_error"
    )
  }

  standard_errors <- as.double(survey::SE(estimates))
  named_values <- stats::setNames(standard_errors^2, result_labels)[labels]
  domain_estimates <- stats::setNames(as.double(estimates$y), result_labels)[labels]
  if (length(named_values) != length(labels) ||
      any(!is.finite(named_values)) || any(named_values < 0)) {
    .bhf_vhat_abort(
      "Taylor sampling variances must be finite and nonnegative.",
      "bhf_vhat_value_error"
    )
  }

  list(
    enabled = TRUE,
    a_star_available = TRUE,
    stan_values = unname(named_values),
    named_values = named_values,
    provenance = list(
      mode = "taylor",
      source = "builtin_survey_taylor",
      supplied_method = NULL,
      engine = "survey",
      engine_version = as.character(utils::packageVersion("survey")),
      estimator = "Hajek domain proportion",
      variance_method = "Taylor linearization",
      design = "stratified one-stage PSU with-replacement approximation",
      weights = "raw",
      fpc = FALSE,
      lonely_psu = "fail",
      adjust_domain_lonely = FALSE,
      options_restoration_registered = TRUE,
      n = nrow(survey_data),
      n_domains = length(labels),
      n_strata = length(unique(survey_data$stratum_id)),
      n_psus = length(unique(survey_data$psu_flat_id)),
      domain_estimates = domain_estimates,
      units = "probability_squared",
      fixed_input = TRUE,
      uncertainty_propagated = FALSE,
      stan_placeholder = FALSE,
      domain_labels = labels
    )
  )
}

#' Resolve the de-attenuation mode
#'
#' Internal compatibility helper. A `NULL` new-style argument means that the
#' caller did not explicitly select a mode, so the article-aligned Taylor mode
#' is the default. The deprecated logical argument is retained for the 0.4.x
#' transition only.
#'
#' @keywords internal
resolve_deattenuation_mode <- function(deattenuation = NULL,
                                       sampling_variances = NULL,
                                       use_deattenuation = NULL) {
  choices <- c("taylor", "supplied", "none")

  if (!is.null(use_deattenuation)) {
    if (!is.null(deattenuation)) {
      .bhf_vhat_abort(
        "Specify only one of deattenuation and deprecated use_deattenuation.",
        "bhf_vhat_mode_error"
      )
    }
    if (!is.logical(use_deattenuation) || length(use_deattenuation) != 1L ||
        is.na(use_deattenuation)) {
      .bhf_vhat_abort(
        "use_deattenuation must be one non-missing logical value.",
        "bhf_vhat_mode_error"
      )
    }

    .bhf_vhat_warn(
      paste0(
        "use_deattenuation is deprecated in bhfvar 0.4.x; use deattenuation ",
        "instead. It may be removed no earlier than 0.5.0, and the old ",
        "heuristic sampling-variance calculation is not retained."
      ),
      "bhf_deprecated_argument_warning"
    )

    if (!use_deattenuation) {
      return("none")
    }
    if (!is.null(sampling_variances)) {
      return("supplied")
    }
    return("taylor")
  }

  if (is.null(deattenuation)) {
    return("taylor")
  }

  tryCatch(
    match.arg(deattenuation, choices),
    error = function(error) {
      .bhf_vhat_abort(
        "deattenuation must be exactly one of: taylor, supplied, none.",
        "bhf_vhat_mode_error"
      )
    }
  )
}

#' Resolve design-based sampling variances
#'
#' Validates and aligns fixed, externally calculated domain sampling variances,
#' or records that de-attenuation is unavailable. Taylor calculation is
#' delegated to `estimate_taylor_vhat()`.
#'
#' @keywords internal
resolve_sampling_variances <- function(deattenuation = NULL,
                                       sampling_variances = NULL,
                                       sampling_variance_method = NULL,
                                       domain_labels,
                                       analysis_data = NULL,
                                       use_deattenuation = NULL) {
  labels <- .validate_domain_labels(domain_labels)
  mode <- resolve_deattenuation_mode(
    deattenuation = deattenuation,
    sampling_variances = sampling_variances,
    use_deattenuation = use_deattenuation
  )

  if (identical(mode, "none")) {
    if (!is.null(sampling_variances) || !is.null(sampling_variance_method)) {
      .bhf_vhat_abort(
        paste0(
          "sampling_variances and sampling_variance_method must be NULL ",
          "when deattenuation = 'none'."
        ),
        "bhf_vhat_mode_error"
      )
    }

    return(list(
      enabled = FALSE,
      a_star_available = FALSE,
      stan_values = rep(0, length(labels)),
      named_values = NULL,
      provenance = list(
        mode = "none",
        source = "none",
        supplied_method = NULL,
        units = "probability_squared",
        fixed_input = FALSE,
        uncertainty_propagated = FALSE,
        stan_placeholder = TRUE,
        domain_labels = labels
      )
    ))
  }

  if (identical(mode, "supplied")) {
    aligned <- .validate_supplied_vhat(sampling_variances, labels)
    supplied_method <- .validate_supplied_vhat_method(
      sampling_variance_method
    )

    return(list(
      enabled = TRUE,
      a_star_available = TRUE,
      stan_values = unname(aligned),
      named_values = aligned,
      provenance = list(
        mode = "supplied",
        source = "external",
        supplied_method = supplied_method,
        units = "probability_squared",
        fixed_input = TRUE,
        uncertainty_propagated = FALSE,
        stan_placeholder = FALSE,
        domain_labels = labels
      )
    ))
  }

  if (!is.null(sampling_variances) || !is.null(sampling_variance_method)) {
    .bhf_vhat_abort(
      paste0(
        "sampling_variances and sampling_variance_method must be NULL ",
        "when deattenuation = 'taylor'."
      ),
      "bhf_vhat_mode_error"
    )
  }
  estimate_taylor_vhat(
    analysis_data = analysis_data,
    domain_labels = labels
  )
}
