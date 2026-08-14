.bhf_legacy_fixture <- function(name) {
  dget(testthat::test_path("..", "fixtures", name))
}

.bhf_current_data_fixture <- function() {
  fixture <- make_tiny_crossed_design_fixture()
  suppressWarnings(prepare_bhf_data(
    fixture$data,
    outcome = "outcome",
    domain = "state",
    strata = "stratum",
    psu = "psu",
    weights = "weight",
    population_shares = fixture$known$population_shares,
    deattenuation = "supplied",
    sampling_variances = fixture$known$vhat,
    sampling_variance_method = "external_taylor"
  ))
}

test_that("text fixtures freeze deterministic 0.3.0 data and fit shapes", {
  data_path <- testthat::test_path(
    "..", "fixtures", "bhf_data-0.3.0.dput"
  )
  fit_path <- testthat::test_path(
    "..", "fixtures", "bhf_fit-0.3.0.dput"
  )
  legacy_data <- dget(data_path)
  legacy_fit <- dget(fit_path)

  expect_s3_class(legacy_data, "bhf_data")
  expect_s3_class(legacy_fit, "bhf_fit")
  expect_identical(legacy_fit$data, legacy_data)
  expect_false(any(c(
    "schema_version", "contract_id", "stan_data_schema_version"
  ) %in% names(legacy_data)))
  expect_false(any(c(
    "data_schema_version", "psu_flat_id", "sigma_state_prior_code"
  ) %in% names(legacy_data$stan_data)))
  expect_identical(dget(data_path), legacy_data)
  expect_identical(dget(fit_path), legacy_fit)
})

test_that("schema detection distinguishes current, legacy, and unsupported", {
  legacy_data <- .bhf_legacy_fixture("bhf_data-0.3.0.dput")
  legacy_fit <- .bhf_legacy_fixture("bhf_fit-0.3.0.dput")
  current_data <- .bhf_current_data_fixture()
  current_fit <- structure(
    list(
      schema_version = "0.5.0",
      contract_id = "bhfvar-fit-contract-0.5.0",
      data = current_data
    ),
    class = c("bhf_fit", "list")
  )

  data_detection <- detect_bhf_object_schema(legacy_data)
  fit_detection <- detect_bhf_object_schema(legacy_fit)
  expect_s3_class(data_detection, "bhf_schema_detection")
  expect_identical(
    unname(unlist(data_detection[c(
      "object_type", "schema_version", "status", "inferred"
    )])),
    c("bhf_data", "0.3.0", "legacy", "TRUE")
  )
  expect_identical(fit_detection$object_type, "bhf_fit")
  expect_identical(fit_detection$schema_version, "0.3.0")
  expect_identical(fit_detection$status, "legacy")
  expect_true(fit_detection$inferred)

  expect_identical(detect_bhf_object_schema(current_data)$status, "current")
  expect_identical(detect_bhf_object_schema(current_fit)$status, "current")
  unversioned_fit <- structure(
    list(data = current_data), class = c("bhf_fit", "list")
  )
  expect_identical(
    detect_bhf_object_schema(unversioned_fit)$status, "unsupported"
  )
  expect_identical(detect_bhf_object_schema(list())$status, "not_bhf")

  unsupported <- current_data
  unsupported$schema_version <- "9.9.9"
  unsupported_detection <- detect_bhf_object_schema(unsupported)
  expect_identical(unsupported_detection$status, "unsupported")
  expect_identical(unsupported_detection$schema_version, "9.9.9")
})

test_that("legacy objects fail closed with exact re-prepare and refit actions", {
  legacy_data <- .bhf_legacy_fixture("bhf_data-0.3.0.dput")
  legacy_fit <- .bhf_legacy_fixture("bhf_fit-0.3.0.dput")
  data_before <- legacy_data
  fit_before <- legacy_fit

  data_error <- tryCatch(
    .bhf_legacy_assert_current_bhf_data(legacy_data),
    error = identity
  )
  fit_error <- tryCatch(
    .bhf_legacy_assert_current_bhf_fit(legacy_fit),
    error = identity
  )

  expect_s3_class(data_error, "bhf_legacy_object_error")
  expect_s3_class(data_error, "bhf_schema_error")
  expect_identical(data_error$detected_schema, "0.3.0")
  expect_identical(data_error$current_schema, "0.5.0")
  expect_identical(data_error$required_action, "reprepare")
  expect_false(data_error$automatic_adapter)
  expect_match(conditionMessage(data_error), "model and estimand semantics changed")
  expect_match(conditionMessage(data_error), "prepare_bhf_data\\(\\)")
  expect_match(conditionMessage(data_error), "silent reinterpretation")

  expect_s3_class(fit_error, "bhf_legacy_object_error")
  expect_identical(fit_error$required_action, "refit")
  expect_match(conditionMessage(fit_error), "bhf_fit\\(\\)")
  expect_identical(legacy_data, data_before)
  expect_identical(legacy_fit, fit_before)
})

test_that("current assertions pass and unknown schemas fail closed", {
  current_data <- .bhf_current_data_fixture()
  current_fit <- structure(
    list(
      schema_version = "0.5.0",
      contract_id = "bhfvar-fit-contract-0.5.0",
      data = current_data
    ),
    class = c("bhf_fit", "list")
  )
  expect_invisible(.bhf_legacy_assert_current_bhf_data(current_data))
  expect_invisible(.bhf_legacy_assert_current_bhf_fit(current_fit))

  unknown <- current_data
  unknown$schema_version <- "9.9.9"
  error <- tryCatch(
    .bhf_legacy_assert_current_bhf_data(unknown), error = identity
  )
  expect_s3_class(error, "bhf_unsupported_schema_error")
  expect_identical(error$detected_schema, "9.9.9")
  expect_identical(error$required_action, "reprepare")
})

test_that("old-name inventory is explicit and is never an adaptation map", {
  inventory <- .bhf_legacy_old_name_inventory()
  legacy_data <- .bhf_legacy_fixture("bhf_data-0.3.0.dput")
  legacy_fit <- .bhf_legacy_fixture("bhf_fit-0.3.0.dput")

  expect_setequal(names(legacy_data), inventory$data_top_level)
  expect_setequal(names(legacy_data$stan_data), inventory$stan_data)
  expect_setequal(names(legacy_fit), inventory$fit_top_level)
  expect_identical(
    legacy_fit$stanfit$draw_names,
    inventory$posterior_draws
  )
  expect_true(all(c(
    "icc_state", "p_state_marginal", "icc_prob", "icc_deatten",
    "reliability_state", "log_lik"
  ) %in% inventory$posterior_draws))
  expect_false(any(c(
    "icc_state_latent", "p_state_A", "prop_between_B",
    "prop_between_A_star", "log_lik_raw", "log_lik_pseudo"
  ) %in% inventory$posterior_draws))
})

test_that("legacy detection has a one-release fail-closed policy", {
  policy <- .bhf_legacy_policy()
  expect_identical(policy$legacy_schema_version, "0.3.0")
  expect_identical(policy$current_schema_version, "0.5.0")
  expect_identical(policy$detection_support_release, "0.5.0")
  expect_identical(policy$removal_eligible_release, "0.5.0")
  expect_identical(policy$support_window, "one_minor_release")
  expect_identical(policy$behavior, "detect_and_fail_closed")
  expect_false(policy$automatic_adapter)
  expect_false(policy$silent_reinterpretation)
  expect_identical(policy$data_action, "reprepare")
  expect_identical(policy$fit_action, "refit")
})
