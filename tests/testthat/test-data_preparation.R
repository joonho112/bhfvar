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

  # Complete inputs retain a reversible identity row mapping.
  expect_identical(
    prepared$input_info$retained_rows,
    seq_len(nrow(bhf_synthetic_data))
  )
  expect_identical(prepared$input_info$dropped_rows, integer())
  expect_equal(nrow(prepared$input_info$missing_reason_ledger), 0)
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

test_that("prepare_bhf_data uses one complete-case analysis universe", {
  n <- 60
  df <- data.frame(
    outcome = rep(c(0, 1), n / 2),
    state = rep(c("A", "B", "C"), each = n / 3),
    stratum = rep(seq_len(6), each = n / 6),
    psu = seq_len(n),
    weight = as.double(seq_len(n))
  )

  df$outcome[2] <- NA
  df$state[17] <- NA
  df$stratum[33] <- NA
  df$psu[45] <- NA
  df$weight[c(2, 59)] <- NA
  dropped_rows <- c(2L, 17L, 33L, 45L, 59L)

  expect_warning(
    prepared <- prepare_bhf_data(
      df, "outcome", "state", "stratum", "psu", "weight"
    ),
    paste0(
      "Excluded 5 observations with missing required analysis values ",
      "\\(outcome=1, domain=1, strata=1, psu=1, weights=2\\)"
    )
  )

  expect_equal(prepared$input_info$n_original, n)
  expect_equal(prepared$input_info$n_used, n - length(dropped_rows))
  expect_identical(prepared$input_info$dropped_rows, dropped_rows)
  expect_identical(
    prepared$input_info$retained_rows,
    setdiff(seq_len(n), dropped_rows)
  )
  expect_equal(prepared$stan_data$N, n - length(dropped_rows))
  expect_false(anyNA(prepared$stan_data$y))
  expect_false(anyNA(prepared$stan_data$state_id))
  expect_false(anyNA(prepared$stan_data$stratum_id))
  expect_false(anyNA(prepared$stan_data$psu_in_stratum_id))
  expect_false(anyNA(prepared$stan_data$w_lik))

  expected_ledger <- data.frame(
    original_row = c(2L, 2L, 17L, 33L, 45L, 59L),
    original_row_name = as.character(c(2L, 2L, 17L, 33L, 45L, 59L)),
    field = c("outcome", "weights", "domain", "strata", "psu", "weights"),
    variable = c("outcome", "weight", "state", "stratum", "psu", "weight"),
    stringsAsFactors = FALSE
  )
  expect_identical(prepared$input_info$missing_reason_ledger, expected_ledger)
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
