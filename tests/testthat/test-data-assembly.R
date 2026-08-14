test_that("crossed supplied-vhat data assembly round-trips the frozen contract", {
  fixture <- make_tiny_crossed_design_fixture()
  supplied_vhat <- fixture$known$vhat[c("C", "A", "B")]

  expect_warning(
    prepared <- prepare_bhf_data(
      fixture$data,
      outcome = "outcome",
      domain = "state",
      strata = "stratum",
      psu = "psu",
      weights = "weight",
      population_shares = fixture$known$population_shares,
      deattenuation = "supplied",
      sampling_variances = supplied_vhat,
      sampling_variance_method = "external_taylor"
    ),
    "Sample size \\(12\\) is very small"
  )

  expect_identical(prepared$schema_version, "0.5.0")
  expect_identical(prepared$contract_id, "bhfvar-data-contract-0.5.0")
  expect_identical(prepared$stan_data_schema_version, 1L)
  expect_identical(prepared$stan_data$data_schema_version, 1L)
  expect_identical(prepared$input_info$schema_version, "0.5.0")
  expect_identical(prepared$input_info$contract_id,
                   "bhfvar-data-contract-0.5.0")

  truth <- fixture$truth
  stan_data <- prepared$stan_data
  expect_identical(
    c(N = stan_data$N, S = stan_data$S, H = stan_data$H, J = stan_data$J),
    truth$dimensions
  )
  expect_identical(as.vector(stan_data$y), fixture$data$outcome)
  expect_identical(as.vector(stan_data$state_id), unname(truth$state_id))
  expect_identical(as.vector(stan_data$stratum_id), unname(truth$stratum_id))
  expect_identical(as.vector(stan_data$J_h), unname(truth$J_h))
  expect_identical(as.vector(stan_data$psu_start), unname(truth$psu_start))
  expect_identical(
    as.vector(stan_data$psu_in_stratum_id),
    unname(truth$psu_within_stratum_id)
  )
  expect_identical(as.vector(stan_data$psu_flat_id),
                   unname(truth$psu_flat_id))

  expect_identical(prepared$mapping$domain$label, truth$state_levels)
  expect_identical(prepared$mapping$domain$id, 1:3)
  expect_identical(prepared$mapping$stratum$label, truth$stratum_levels)
  expect_identical(prepared$mapping$stratum$id, 1:2)
  expect_identical(prepared$mapping$psu$label, c("P1", "P2", "P1", "P2"))
  expect_identical(prepared$mapping$psu$id, 1:4)
  expect_identical(prepared$mapping$psu$stratum_id, c(1L, 1L, 2L, 2L))

  expect_identical(prepared$row_provenance$n_original, 12L)
  expect_identical(prepared$row_provenance$n_used, 12L)
  expect_identical(prepared$row_provenance$retained_rows, 1:12)
  expect_identical(prepared$row_provenance$dropped_rows, integer())
  expect_equal(nrow(prepared$row_provenance$missing_reason_ledger), 0)
  expect_identical(prepared$provenance$rows, prepared$row_provenance)

  analysis_data <- prepared$analysis_data
  expect_identical(analysis_data$original_row, 1:12)
  expect_identical(analysis_data$y, fixture$data$outcome)
  expect_identical(analysis_data$domain_label, fixture$data$state)
  expect_identical(analysis_data$state_id, unname(truth$state_id))
  expect_identical(analysis_data$stratum_label, fixture$data$stratum)
  expect_identical(analysis_data$stratum_id, unname(truth$stratum_id))
  expect_identical(analysis_data$psu_label, fixture$data$psu)
  expect_identical(analysis_data$psu_in_stratum_id,
                   unname(truth$psu_within_stratum_id))
  expect_identical(analysis_data$psu_flat_id, unname(truth$psu_flat_id))
  expect_identical(analysis_data$raw_weight, fixture$data$weight)
  expect_equal(analysis_data$w_lik, unname(truth$mean_one_weights),
               tolerance = 1e-12)

  expect_identical(prepared$weight_info$original_row, 1:12)
  expect_identical(prepared$weight_info$raw, fixture$data$weight)
  expect_equal(prepared$weight_info$likelihood,
               unname(truth$mean_one_weights), tolerance = 1e-12)
  expect_identical(prepared$provenance$weights,
                   prepared$weight_info$provenance)
  expect_identical(prepared$provenance$weights$method, "mean_one")
  expect_equal(prepared$provenance$weights$raw_sum, 78)
  expect_equal(prepared$provenance$weights$likelihood_sum, 12,
               tolerance = 1e-12)

  expect_identical(prepared$population_share_info$values,
                   fixture$known$population_shares)
  expect_identical(as.vector(stan_data$w_state_pop_share),
                   unname(fixture$known$population_shares))
  expect_identical(prepared$provenance$population_shares$source,
                   "external_known")
  expect_identical(prepared$provenance$population_shares$domain_labels,
                   truth$state_levels)

  expect_identical(prepared$sampling_variance_info$named_values,
                   fixture$known$vhat)
  expect_identical(as.vector(stan_data$vhat_state),
                   unname(fixture$known$vhat))
  expect_identical(prepared$provenance$sampling_variances,
                   prepared$sampling_variance_info$provenance)
  expect_identical(prepared$provenance$sampling_variances$mode, "supplied")
  expect_identical(
    prepared$provenance$sampling_variances$supplied_method,
    "external_taylor"
  )

  expected_alpha_mean <- stats::qlogis(
    sum(fixture$data$weight * fixture$data$outcome) /
      sum(fixture$data$weight)
  )
  expect_equal(stan_data$prior_alpha_mean, expected_alpha_mean,
               tolerance = 1e-12)
  expect_identical(stan_data$prior_alpha_sd, 0.5)
  expect_identical(prepared$prior_info$alpha$mean_source,
                   "estimated_from_raw_weighted_outcome_logit")
  expect_identical(prepared$prior_info$alpha$sd_source, "default")
  expect_equal(
    prepared$prior_info$alpha$design_weighted_prevalence,
    sum(fixture$data$weight * fixture$data$outcome) /
      sum(fixture$data$weight),
    tolerance = 1e-12
  )
  expect_identical(
    prepared$prior_info$random_effect_sd$family,
    "half_student_t"
  )
  expect_identical(prepared$provenance$prior, prepared$prior_info)
})

test_that("data assembly preserves transition fields alongside new surfaces", {
  fixture <- make_tiny_crossed_design_fixture()
  prepared <- suppressWarnings(prepare_bhf_data(
    fixture$data,
    "outcome", "state", "stratum", "psu", "weight",
    population_shares = fixture$known$population_shares,
    deattenuation = "supplied",
    sampling_variances = fixture$known$vhat,
    sampling_variance_method = "external_other",
    prior_alpha_mean = -1.25,
    prior_alpha_sd = 0.75
  ))

  expect_true(all(c(
    "stan_data", "mapping", "weight_info", "population_share_info",
    "sampling_variance_info", "analysis_data", "domain_summary", "input_info"
  ) %in% names(prepared)))
  expect_identical(prepared$input_info$retained_rows,
                   prepared$row_provenance$retained_rows)
  expect_identical(prepared$input_info$dropped_rows,
                   prepared$row_provenance$dropped_rows)
  expect_identical(prepared$prior_info$alpha$mean, -1.25)
  expect_identical(prepared$prior_info$alpha$mean_source, "user_supplied")
  expect_identical(prepared$prior_info$alpha$sd, 0.75)
  expect_identical(prepared$prior_info$alpha$sd_source, "user_supplied")
})

test_that("automatic intercept prior fails closed for degenerate outcomes", {
  fixture <- make_tiny_crossed_design_fixture()
  fixture$data$outcome <- 0L

  expect_error(
    suppressWarnings(prepare_bhf_data(
      fixture$data,
      "outcome", "state", "stratum", "psu", "weight",
      population_shares = fixture$known$population_shares,
      deattenuation = "none"
    )),
    "strictly between 0 and 1"
  )
  expect_s3_class(
    suppressWarnings(prepare_bhf_data(
      fixture$data,
      "outcome", "state", "stratum", "psu", "weight",
      population_shares = fixture$known$population_shares,
      deattenuation = "none",
      prior_alpha_mean = -4
    )),
    "bhf_data"
  )
})
