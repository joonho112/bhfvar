test_that("mean-one likelihood weights have exact expected values and ratios", {
  raw <- c(a = 1, b = 2, c = 3, d = 4)

  scaled <- scale_likelihood_weights(raw)

  expect_equal(scaled, c(a = 0.4, b = 0.8, c = 1.2, d = 1.6),
               tolerance = 1e-12)
  expect_equal(sum(scaled), length(raw), tolerance = 1e-12)
  expect_equal(mean(scaled), 1, tolerance = 1e-12)
  expect_equal(unname(scaled / scaled[[1]]), unname(raw / raw[[1]]),
               tolerance = 1e-12)
})

test_that("mean-one scaling is scale- and permutation-invariant", {
  raw <- c(0.2, 3.7, 10.5, 1.1, 8.3)
  permutation <- c(5L, 2L, 4L, 1L, 3L)
  reference <- scale_likelihood_weights(raw)

  expect_equal(scale_likelihood_weights(raw * 1e12), reference,
               tolerance = 1e-12)
  expect_equal(scale_likelihood_weights(raw * 1e-12), reference,
               tolerance = 1e-12)
  expect_equal(
    scale_likelihood_weights(raw[permutation])[order(permutation)],
    reference,
    tolerance = 1e-12
  )
})

test_that("likelihood weight scaling rejects invalid raw inputs", {
  expect_error(scale_likelihood_weights(character()), "non-empty numeric")
  expect_error(scale_likelihood_weights(numeric()), "non-empty numeric")
  expect_error(scale_likelihood_weights(c(1, NA_real_)), "finite, strictly positive")
  expect_error(scale_likelihood_weights(c(1, Inf)), "finite, strictly positive")
  expect_error(scale_likelihood_weights(c(1, 0)), "finite, strictly positive")
  expect_error(scale_likelihood_weights(c(1, -1)), "finite, strictly positive")
  expect_error(scale_likelihood_weights(c(1, 2), method = "unknown"),
               "should be one of")
})

test_that("legacy D2 is explicit, warned, and globally normalized", {
  raw <- c(1, 3, 2, 6)
  state_id <- c(1L, 1L, 2L, 2L)

  expect_warning(
    scaled <- scale_likelihood_weights(raw, "legacy_d2", state_id),
    "legacy_d2.*deprecated",
    class = "bhf_legacy_scaling_warning"
  )

  expect_equal(scaled, c(0.5, 1.5, 0.5, 1.5), tolerance = 1e-12)
  expect_equal(sum(scaled), length(raw), tolerance = 1e-12)
  expect_equal(scaled[2] / scaled[1], raw[2] / raw[1], tolerance = 1e-12)
  expect_equal(scaled[4] / scaled[3], raw[4] / raw[3], tolerance = 1e-12)
})

test_that("legacy D2 validates domains and remains invariant", {
  raw <- c(1, 3, 2, 6, 5, 4)
  state_id <- c("A", "A", "B", "B", "C", "C")
  permutation <- c(6L, 1L, 4L, 2L, 5L, 3L)

  expect_warning(
    reference <- scale_likelihood_weights(raw, "legacy_d2", state_id),
    "legacy_d2.*deprecated"
  )
  expect_warning(
    rescaled <- scale_likelihood_weights(raw * 100, "legacy_d2", state_id),
    "legacy_d2.*deprecated"
  )
  expect_warning(
    permuted <- scale_likelihood_weights(
      raw[permutation], "legacy_d2", state_id[permutation]
    ),
    "legacy_d2.*deprecated"
  )

  expect_equal(rescaled, reference, tolerance = 1e-12)
  expect_equal(permuted[order(permutation)], reference, tolerance = 1e-12)

  expect_warning(
    expect_error(
      scale_likelihood_weights(raw, "legacy_d2"),
      "state_id must be an aligned"
    ),
    "legacy_d2.*deprecated"
  )
  expect_warning(
    expect_error(
      scale_likelihood_weights(raw, "legacy_d2", state_id[-1]),
      "state_id must be an aligned"
    ),
    "legacy_d2.*deprecated"
  )
  expect_warning(
    expect_error(
      scale_likelihood_weights(raw, "legacy_d2", replace(state_id, 1, NA)),
      "state_id must be an aligned"
    ),
    "legacy_d2.*deprecated"
  )
})

test_that("prepare_bhf_data defaults to canonical mean-one weights", {
  fixture <- make_tiny_crossed_design_fixture()

  expect_warning(
    prepared <- prepare_bhf_data(
      fixture$data,
      outcome = "outcome",
      domain = "state",
      strata = "stratum",
      psu = "psu",
      weights = "weight",
      population_shares = fixture$known$population_shares
    ),
    "Sample size \\(12\\) is very small"
  )

  expect_equal(as.vector(prepared$stan_data$w_lik),
               unname(fixture$truth$mean_one_weights), tolerance = 1e-12)
  expect_equal(sum(prepared$stan_data$w_lik), prepared$stan_data$N,
               tolerance = 1e-12)
  expect_identical(prepared$input_info$weight_scaling, "mean_one")
  expect_identical(prepared$weight_info$original_row, seq_len(12))
  expect_identical(prepared$weight_info$raw, fixture$data$weight)
  expect_equal(prepared$weight_info$likelihood,
               unname(fixture$truth$mean_one_weights), tolerance = 1e-12)
  expect_identical(prepared$weight_info$provenance$method, "mean_one")
  expect_identical(
    prepared$weight_info$provenance$normalization,
    "global_sum_N_in_R"
  )
  expect_equal(prepared$weight_info$provenance$raw_sum,
               sum(fixture$data$weight))
  expect_equal(prepared$weight_info$provenance$likelihood_sum, 12,
               tolerance = 1e-12)
  expect_equal(prepared$weight_info$provenance$global_factor, 2 / 13,
               tolerance = 1e-12)
  expect_null(prepared$weight_info$provenance$legacy_domain_factors)
})

test_that("prepare_bhf_data exposes warned legacy D2 provenance", {
  fixture <- make_tiny_crossed_design_fixture()

  expect_warning(
    expect_warning(
      prepared <- prepare_bhf_data(
        fixture$data,
        outcome = "outcome",
        domain = "state",
        strata = "stratum",
        psu = "psu",
        weights = "weight",
        population_shares = fixture$known$population_shares,
        weight_scaling = "legacy_d2"
      ),
      "legacy_d2.*deprecated"
    ),
    "Sample size \\(12\\) is very small"
  )

  expect_identical(prepared$input_info$weight_scaling, "legacy_d2")
  expect_identical(prepared$weight_info$provenance$method, "legacy_d2")
  expect_equal(sum(prepared$weight_info$likelihood), 12, tolerance = 1e-12)
  expect_named(
    prepared$weight_info$provenance$legacy_domain_factors,
    fixture$truth$state_levels
  )
  expect_false(isTRUE(all.equal(
    prepared$weight_info$likelihood,
    unname(fixture$truth$mean_one_weights),
    tolerance = 1e-12
  )))
})
