# ==============================================================================
# Tests for result extractor functions
# ==============================================================================
#
# These tests verify the output structure of extractor functions.
# Full integration tests require a fitted model.
#
# Author: JoonHo Lee
# Date: 2025-01-24
# ==============================================================================

test_that("variance_decomposition requires bhf_fit object", {
  expect_error(
    variance_decomposition(list(x = 1)),
    "bhf_fit"
  )
})

test_that("domain_estimates requires bhf_fit object", {
  expect_error(
    domain_estimates(list(x = 1)),
    "bhf_fit"
  )
})

test_that("overall_estimate requires bhf_fit object", {
  expect_error(
    overall_estimate(list(x = 1)),
    "bhf_fit"
  )
})

test_that("log_lik requires bhf_fit object", {
  expect_error(
    log_lik(list(x = 1)),
    "bhf_fit"
  )
})

test_that("domain_estimates accepts valid type argument", {
  expect_error(
    domain_estimates(list(), type = "invalid"),
    "arg"
  )
})
