test_that("de-attenuation mode resolution has explicit defaults", {
  expect_identical(resolve_deattenuation_mode(), "taylor")
  expect_identical(resolve_deattenuation_mode("taylor"), "taylor")
  expect_identical(resolve_deattenuation_mode("supplied"), "supplied")
  expect_identical(resolve_deattenuation_mode("none"), "none")

  expect_error(
    resolve_deattenuation_mode("invalid"),
    class = "bhf_vhat_mode_error"
  )
  expect_error(
    resolve_deattenuation_mode(c("supplied", "none")),
    class = "bhf_vhat_mode_error"
  )
})

test_that("supplied sampling variances round-trip in canonical domain order", {
  result <- resolve_sampling_variances(
    deattenuation = "supplied",
    sampling_variances = c(B = 0.04, A = 0.01),
    sampling_variance_method = "external_replicate",
    domain_labels = c("A", "B")
  )

  expect_true(result$enabled)
  expect_true(result$a_star_available)
  expect_identical(result$named_values, c(A = 0.01, B = 0.04))
  expect_identical(result$stan_values, c(0.01, 0.04))
  expect_identical(result$provenance$mode, "supplied")
  expect_identical(result$provenance$source, "external")
  expect_identical(
    result$provenance$supplied_method,
    "external_replicate"
  )
  expect_identical(result$provenance$units, "probability_squared")
  expect_true(result$provenance$fixed_input)
  expect_false(result$provenance$uncertainty_propagated)
  expect_false(result$provenance$stan_placeholder)
})

test_that("permuted supplied names are aligned and zero is valid", {
  first <- resolve_sampling_variances(
    deattenuation = "supplied",
    sampling_variances = c(A = 0, B = 0.02),
    sampling_variance_method = "external_taylor",
    domain_labels = c("A", "B")
  )
  permuted <- resolve_sampling_variances(
    deattenuation = "supplied",
    sampling_variances = c(B = 0.02, A = 0),
    sampling_variance_method = "external_taylor",
    domain_labels = c("A", "B")
  )

  expect_identical(first$named_values, c(A = 0, B = 0.02))
  expect_identical(permuted$named_values, first$named_values)
  expect_identical(permuted$stan_values, c(0, 0.02))
})

test_that("supplied sampling variances reject invalid values", {
  invalid_values <- list(
    c(A = -0.01, B = 0.02),
    c(A = NA_real_, B = 0.02),
    c(A = NaN, B = 0.02),
    c(A = Inf, B = 0.02)
  )

  for (values in invalid_values) {
    expect_error(
      resolve_sampling_variances(
        deattenuation = "supplied",
        sampling_variances = values,
        sampling_variance_method = "external_other",
        domain_labels = c("A", "B")
      ),
      class = "bhf_vhat_value_error"
    )
  }

  expect_error(
    resolve_sampling_variances(
      deattenuation = "supplied",
      sampling_variances = c(A = "0.01", B = "0.02"),
      sampling_variance_method = "external_other",
      domain_labels = c("A", "B")
    ),
    class = "bhf_vhat_value_error"
  )
})

test_that("supplied sampling variances require exact names", {
  invalid_names <- list(
    c(0.01, 0.02),
    stats::setNames(c(0.01, 0.02), c("", "B")),
    stats::setNames(c(0.01, 0.02), c("A", "A")),
    c(A = 0.01),
    c(A = 0.01, B = 0.02, C = 0.03),
    c(A = 0.01, C = 0.02)
  )

  for (values in invalid_names) {
    expect_error(
      resolve_sampling_variances(
        deattenuation = "supplied",
        sampling_variances = values,
        sampling_variance_method = "external_other",
        domain_labels = c("A", "B")
      ),
      class = "bhf_vhat_alignment_error"
    )
  }
})

test_that("supplied mode requires an external method declaration", {
  expect_error(
    resolve_sampling_variances(
      deattenuation = "supplied",
      sampling_variances = c(A = 0.01, B = 0.02),
      domain_labels = c("A", "B")
    ),
    class = "bhf_vhat_mode_error"
  )
  expect_error(
    resolve_sampling_variances(
      deattenuation = "supplied",
      sampling_variances = c(A = 0.01, B = 0.02),
      sampling_variance_method = "kish",
      domain_labels = c("A", "B")
    ),
    class = "bhf_vhat_mode_error"
  )
})

test_that("none mode exposes only a disabled Stan placeholder", {
  result <- resolve_sampling_variances(
    deattenuation = "none",
    domain_labels = c("A", "B", "C")
  )

  expect_false(result$enabled)
  expect_false(result$a_star_available)
  expect_identical(result$stan_values, c(0, 0, 0))
  expect_null(result$named_values)
  expect_identical(result$provenance$mode, "none")
  expect_identical(result$provenance$source, "none")
  expect_true(result$provenance$stan_placeholder)
  expect_false(result$provenance$fixed_input)
  expect_false(result$provenance$uncertainty_propagated)
})

test_that("mode-specific arguments fail closed", {
  expect_error(
    resolve_sampling_variances(
      deattenuation = "supplied",
      domain_labels = c("A", "B"),
      sampling_variance_method = "external_other"
    ),
    class = "bhf_vhat_mode_error"
  )
  expect_error(
    resolve_sampling_variances(
      deattenuation = "none",
      sampling_variances = c(A = 0.01, B = 0.02),
      sampling_variance_method = "external_other",
      domain_labels = c("A", "B")
    ),
    class = "bhf_vhat_mode_error"
  )
  expect_error(
    resolve_sampling_variances(
      deattenuation = "taylor",
      sampling_variances = c(A = 0.01, B = 0.02),
      sampling_variance_method = "external_taylor",
      domain_labels = c("A", "B")
    ),
    class = "bhf_vhat_mode_error"
  )
})

test_that("deprecated use_deattenuation maps with warnings", {
  expect_warning(
    expect_identical(
      resolve_deattenuation_mode(use_deattenuation = FALSE),
      "none"
    ),
    class = "bhf_deprecated_argument_warning"
  )
  expect_warning(
    expect_identical(
      resolve_deattenuation_mode(use_deattenuation = TRUE),
      "taylor"
    ),
    class = "bhf_deprecated_argument_warning"
  )
  expect_warning(
    expect_identical(
      resolve_deattenuation_mode(
        sampling_variances = c(A = 0.01),
        use_deattenuation = TRUE
      ),
      "supplied"
    ),
    class = "bhf_deprecated_argument_warning"
  )

  expect_warning(
    none_result <- resolve_sampling_variances(
      domain_labels = c("A", "B"),
      use_deattenuation = FALSE
    ),
    class = "bhf_deprecated_argument_warning"
  )
  expect_false(none_result$a_star_available)

  expect_warning(
    supplied_result <- resolve_sampling_variances(
      sampling_variances = c(B = 0.02, A = 0.01),
      sampling_variance_method = "external_taylor",
      domain_labels = c("A", "B"),
      use_deattenuation = TRUE
    ),
    class = "bhf_deprecated_argument_warning"
  )
  expect_identical(supplied_result$named_values, c(A = 0.01, B = 0.02))
})

test_that("deprecated and new mode arguments cannot be combined", {
  expect_error(
    resolve_deattenuation_mode(
      deattenuation = "none",
      use_deattenuation = FALSE
    ),
    class = "bhf_vhat_mode_error"
  )
  expect_error(
    resolve_deattenuation_mode(use_deattenuation = NA),
    class = "bhf_vhat_mode_error"
  )
  expect_error(
    resolve_deattenuation_mode(use_deattenuation = c(TRUE, FALSE)),
    class = "bhf_vhat_mode_error"
  )
})

test_that("domain labels are validated before variance alignment", {
  invalid_labels <- list(
    character(),
    c("A", ""),
    c("A", "A"),
    c("A", NA_character_)
  )

  for (labels in invalid_labels) {
    expect_error(
      resolve_sampling_variances(
        deattenuation = "none",
        domain_labels = labels
      ),
      class = "bhf_vhat_alignment_error"
    )
  }
})

test_that("prepare_bhf_data exposes supplied vhat and fixed-input provenance", {
  fixture <- make_tiny_crossed_design_fixture()

  expect_warning(
    prepared <- prepare_bhf_data(
      fixture$data,
      "outcome", "state", "stratum", "psu", "weight",
      population_shares = fixture$known$population_shares,
      deattenuation = "supplied",
      sampling_variances = fixture$known$vhat[c("C", "A", "B")],
      sampling_variance_method = "external_taylor"
    ),
    "Sample size"
  )

  expect_identical(as.vector(prepared$stan_data$vhat_state),
                   unname(fixture$known$vhat))
  expect_identical(prepared$stan_data$use_deattenuation, 1L)
  expect_identical(prepared$sampling_variance_info$named_values,
                   fixture$known$vhat)
  expect_identical(prepared$input_info$deattenuation, "supplied")
  expect_true(prepared$sampling_variance_info$provenance$fixed_input)
  expect_false(
    prepared$sampling_variance_info$provenance$uncertainty_propagated
  )
})

test_that("prepare_bhf_data none mode keeps a transport-only placeholder", {
  fixture <- make_tiny_crossed_design_fixture()
  prepared <- suppressWarnings(prepare_bhf_data(
    fixture$data,
    "outcome", "state", "stratum", "psu", "weight",
    population_shares = fixture$known$population_shares,
    deattenuation = "none"
  ))

  expect_identical(as.vector(prepared$stan_data$vhat_state), c(0, 0, 0))
  expect_identical(prepared$stan_data$use_deattenuation, 0L)
  expect_false(prepared$sampling_variance_info$a_star_available)
  expect_null(prepared$sampling_variance_info$named_values)
  expect_true(prepared$sampling_variance_info$provenance$stan_placeholder)
  expect_true(all(is.na(prepared$domain_summary$vhat)))
})
