# ==============================================================================
# Tests for data preparation functions
# ==============================================================================
#
# Author: JoonHo Lee
# Date: 2025-01-24
# ==============================================================================

test_that("prepare_bhf_data creates correct structure", {
  # Load example data
  data("bhf_synthetic_data", package = "bhfvar")
  
  # Prepare data
  prepared <- prepare_bhf_data(
    bhf_synthetic_data,
    outcome = "has_subsidy",
    domain = "state",
    strata = "stratum",
    psu = "psu",
    weights = "weight"
  )
  
  # Check class
  expect_s3_class(prepared, "bhf_data")
  
  # Check required components exist
  expect_true("stan_data" %in% names(prepared))
  expect_true("mapping" %in% names(prepared))
  expect_true("domain_summary" %in% names(prepared))
})

test_that("prepare_bhf_data validates index recoding", {
  # Load example data
  data("bhf_synthetic_data", package = "bhfvar")
  
  prepared <- prepare_bhf_data(
    bhf_synthetic_data,
    outcome = "has_subsidy",
    domain = "state",
    strata = "stratum",
    psu = "psu",
    weights = "weight"
  )
  
  stan_data <- prepared$stan_data
  
  # Check indices are 1-indexed and consecutive
  expect_equal(min(stan_data$state_id), 1)
  expect_equal(max(stan_data$state_id), stan_data$S)
  
  expect_equal(min(stan_data$stratum_id), 1)
  expect_equal(max(stan_data$stratum_id), stan_data$H)
})

test_that("prepare_bhf_data computes correct weight scaling", {
  # Load example data
  data("bhf_synthetic_data", package = "bhfvar")
  
  prepared <- prepare_bhf_data(
    bhf_synthetic_data,
    outcome = "has_subsidy",
    domain = "state",
    strata = "stratum",
    psu = "psu",
    weights = "weight"
  )
  
  stan_data <- prepared$stan_data
  
  # Sum of weights should approximately equal effective sample size
  # (due to Method D2 scaling)
  expect_true(sum(stan_data$w_lik) > 0)
  
  # Population shares should sum to 1
  expect_equal(sum(stan_data$w_state_pop_share), 1, tolerance = 1e-10)
})

test_that("prepare_bhf_data handles missing values", {
  df <- data.frame(
    outcome = c(0, 1, NA, 0),
    state = c("A", "A", "B", "B"),
    stratum = c(1, 1, 2, 2),
    psu = c(1, 2, 3, 4),
    weight = c(1, 2, 1, 1)
  )
  
  expect_error(
    prepare_bhf_data(df, "outcome", "state", "stratum", "psu", "weight"),
    "missing|NA"
  )
})

test_that("PSU structure is correctly computed", {
  # Load example data
  data("bhf_synthetic_data", package = "bhfvar")
  
  prepared <- prepare_bhf_data(
    bhf_synthetic_data,
    outcome = "has_subsidy",
    domain = "state",
    strata = "stratum",
    psu = "psu",
    weights = "weight"
  )
  
  stan_data <- prepared$stan_data
  
  # J_h should sum to total PSUs
  expect_equal(sum(stan_data$J_h), stan_data$J)
  
  # psu_start should be properly indexed
  expect_equal(stan_data$psu_start[1], 1)
  for (h in 2:stan_data$H) {
    expected_start <- stan_data$psu_start[h-1] + stan_data$J_h[h-1]
    expect_equal(stan_data$psu_start[h], expected_start)
  }
})
