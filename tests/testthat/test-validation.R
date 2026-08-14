# ==============================================================================
# Tests for validation functions
# ==============================================================================
#
# Author: JoonHo Lee
# Date: 2025-01-24
# ==============================================================================

test_that("validate_input_data catches missing variables", {
  df <- data.frame(x = 1:5, y = 1:5)
  
  expect_error(
    validate_input_data(df, "outcome", "state", "stratum", "psu", "weight"),
    "Missing required variables"
  )
})

test_that("validate_input_data catches non-binary outcome", {
  df <- data.frame(
    outcome = c(0, 1, 2, 3),
    state = c("A", "A", "B", "B"),
    stratum = c(1, 1, 2, 2),
    psu = c(1, 2, 3, 4),
    weight = c(1, 1, 1, 1)
  )
  
  expect_error(
    validate_input_data(df, "outcome", "state", "stratum", "psu", "weight"),
    "binary"
  )
})

test_that("validate_input_data catches non-positive weights", {
  df <- data.frame(
    outcome = c(0, 1, 1, 0),
    state = c("A", "A", "B", "B"),
    stratum = c(1, 1, 2, 2),
    psu = c(1, 2, 3, 4),
    weight = c(1, -1, 1, 1)
  )
  
  expect_error(
    validate_input_data(df, "outcome", "state", "stratum", "psu", "weight"),
    "positive"
  )
})

test_that("validate_input_data accepts valid data", {
  n <- 60
  df <- data.frame(
    outcome = rep(c(0, 1), n / 2),
    state = rep(c("A", "B", "C"), each = n / 3),
    stratum = rep(seq_len(6), each = n / 6),
    psu = seq_len(n),
    weight = rep(c(1, 2), n / 2)
  )
  
  expect_silent(
    validate_input_data(df, "outcome", "state", "stratum", "psu", "weight")
  )
})

test_that("calc_eff_n computes correct effective sample size", {
  # Equal weights: eff_n should equal n
  w_equal <- rep(1, 10)
  expect_equal(calc_eff_n(w_equal), 10)
  
  # Unequal weights: eff_n < n
  w_unequal <- c(1, 1, 1, 1, 10)
  eff_n <- calc_eff_n(w_unequal)
  expect_true(eff_n < 5)
  
  # Formula: (sum(w))^2 / sum(w^2)
  expected <- (sum(w_unequal)^2) / sum(w_unequal^2)
  expect_equal(eff_n, expected)
})

test_that("scale_weights_d2 produces correct scaling", {
  # Create simple test data
  weights <- c(1, 2, 3, 4)
  domain_id <- c(1, 1, 2, 2)
  domain_summary <- data.frame(
    state_id = 1:2,
    weight_sum = c(sum(weights[domain_id == 1]), sum(weights[domain_id == 2])),
    eff_n = c(calc_eff_n(weights[domain_id == 1]),
              calc_eff_n(weights[domain_id == 2]))
  )
  
  # Scale weights
  w_scaled <- scale_weights_d2(weights, domain_id, domain_summary)
  
  # Check that sum of scaled weights equals eff_n per domain
  for (d in unique(domain_id)) {
    idx <- domain_id == d
    w_d <- weights[idx]
    w_scaled_d <- w_scaled[idx]
    
    eff_n_d <- calc_eff_n(w_d)
    expect_equal(sum(w_scaled_d), eff_n_d, tolerance = 1e-10)
  }
})
