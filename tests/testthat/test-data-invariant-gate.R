make_contract_test_object <- function(deattenuation = "supplied",
                                      weight_scaling = "mean_one") {
  fixture <- make_tiny_crossed_design_fixture()
  args <- list(
    data = fixture$data,
    outcome = "outcome",
    domain = "state",
    strata = "stratum",
    psu = "psu",
    weights = "weight",
    population_shares = fixture$known$population_shares,
    weight_scaling = weight_scaling,
    deattenuation = deattenuation,
    prior_alpha_mean = -1.5,
    prior_alpha_sd = 0.5
  )
  if (identical(deattenuation, "supplied")) {
    args$sampling_variances <- fixture$known$vhat
    args$sampling_variance_method <- "external_taylor"
  }
  suppressWarnings(do.call(prepare_bhf_data, args))
}

expect_contract_error <- function(object, invariant) {
  error <- tryCatch(
    validate_bhf_data_contract(object),
    error = identity
  )
  expect_s3_class(error, "bhf_data_contract_error")
  expect_identical(error$invariant, invariant)
}

test_that("whole-object contract accepts all supported data branches", {
  supplied <- make_contract_test_object("supplied", "mean_one")
  none <- make_contract_test_object("none", "mean_one")
  legacy <- make_contract_test_object("supplied", "legacy_d2")
  taylor <- make_contract_test_object("taylor", "mean_one")

  expect_invisible(validate_bhf_data_contract(supplied))
  expect_invisible(validate_bhf_data_contract(none))
  expect_invisible(validate_bhf_data_contract(legacy))
  expect_invisible(validate_bhf_data_contract(taylor))
})

test_that("schema, row, and mapping corruptions fail at stable invariants", {
  canonical <- make_contract_test_object()

  broken <- canonical
  broken$schema_version <- "0.3.0"
  expect_contract_error(broken, "schema.version")

  broken <- canonical
  broken$analysis_data <- broken$analysis_data[-1, , drop = FALSE]
  expect_contract_error(broken, "analysis.shape")

  broken <- canonical
  broken$stan_data$y[1] <- 1L - broken$stan_data$y[1]
  expect_contract_error(broken, "outcome.transport")

  broken <- canonical
  broken$row_provenance$retained_rows[1] <-
    broken$row_provenance$retained_rows[2]
  expect_contract_error(broken, "rows.partition")

  broken <- canonical
  bad_ledger <- data.frame(
    original_row = 1L,
    original_row_name = "1",
    field = "outcome",
    variable = "outcome",
    stringsAsFactors = FALSE
  )
  broken$row_provenance$missing_reason_ledger <- bad_ledger
  broken$provenance$rows <- broken$row_provenance
  expect_contract_error(broken, "rows.ledger")

  broken <- canonical
  broken$analysis_data$state_id[1] <- 2L
  expect_contract_error(broken, "mapping.domain")

  broken <- canonical
  broken$mapping$stratum$label[1] <- "BROKEN"
  expect_contract_error(broken, "mapping.stratum")

  broken <- canonical
  broken$mapping$psu$within_stratum_id[1] <- 2L
  expect_contract_error(broken, "mapping.psu")
})

test_that("weight and share corruptions fail before scientific use", {
  canonical <- make_contract_test_object()

  broken <- canonical
  broken$analysis_data$w_lik[1] <- broken$analysis_data$w_lik[1] * 2
  expect_contract_error(broken, "weights.transport")

  broken <- canonical
  broken$weight_info$provenance$raw_sum <-
    broken$weight_info$provenance$raw_sum + 1
  expect_contract_error(broken, "weights.provenance")

  broken <- canonical
  broken$weight_info$likelihood <- rev(broken$weight_info$likelihood)
  broken$analysis_data$w_lik <- broken$weight_info$likelihood
  broken$stan_data$w_lik <- as.array(broken$weight_info$likelihood)
  expect_contract_error(broken, "weights.normalization")

  broken <- canonical
  share_names <- names(broken$population_share_info$values)
  broken$population_share_info$values <- stats::setNames(
    rev(unname(broken$population_share_info$values)),
    share_names
  )
  expect_contract_error(broken, "shares.transport")

  broken <- canonical
  broken$population_share_info$source <- "unknown"
  expect_contract_error(broken, "shares.provenance")
})

test_that("vhat, prior, and summary corruptions fail at stable invariants", {
  canonical <- make_contract_test_object()

  broken <- canonical
  broken$sampling_variance_info$stan_values[1] <-
    broken$sampling_variance_info$stan_values[1] * 2
  expect_contract_error(broken, "vhat.transport")

  broken <- canonical
  broken$sampling_variance_info$enabled <- FALSE
  expect_contract_error(broken, "vhat.mode")

  broken <- canonical
  broken$sampling_variance_info$provenance$source <- "kish"
  expect_contract_error(broken, "vhat.provenance")

  broken <- canonical
  broken$prior_info$alpha$sd <- 0.75
  expect_contract_error(broken, "prior.provenance")

  broken <- canonical
  broken$prior_info$alpha$sd <- 0.75
  broken$provenance$prior <- broken$prior_info
  expect_contract_error(broken, "prior.transport")

  broken <- canonical
  broken$domain_summary$n[1] <- broken$domain_summary$n[1] + 1L
  expect_contract_error(broken, "summary.transport")
})

test_that("none-mode transport placeholder cannot masquerade as vhat", {
  none <- make_contract_test_object("none")

  broken <- none
  broken$stan_data$vhat_state[1] <- 0.01
  expect_contract_error(broken, "vhat.transport")

  broken <- none
  broken$domain_summary$vhat[1] <- 0
  expect_contract_error(broken, "vhat.mode")
})
