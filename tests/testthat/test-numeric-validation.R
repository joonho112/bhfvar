make_valid_validation_data <- function() {
  n <- 60L
  data.frame(
    y = rep(c(0L, 1L), n / 2L),
    domain = rep(c("A", "B", "C"), each = n / 3L),
    stratum = rep(c("H1", "H2"), each = n / 2L),
    psu = sprintf("P%02d", seq_len(n)),
    weight = rep(c(1, 2, 4), length.out = n),
    stringsAsFactors = FALSE
  )
}

make_valid_stan_data <- function() {
  list(
    data_schema_version = 1L,
    N = 6L,
    S = 3L,
    H = 2L,
    J = 4L,
    y = c(0L, 1L, 0L, 1L, 1L, 0L),
    state_id = c(1L, 2L, 3L, 1L, 2L, 3L),
    stratum_id = c(1L, 1L, 1L, 2L, 2L, 2L),
    J_h = c(2L, 2L),
    psu_start = c(1L, 3L),
    psu_in_stratum_id = c(1L, 1L, 2L, 1L, 2L, 2L),
    psu_flat_id = c(1L, 1L, 2L, 3L, 4L, 4L),
    w_lik = rep(1, 6),
    w_state_pop_share = c(0.2, 0.3, 0.5),
    prior_alpha_mean = -1.5,
    prior_alpha_sd = 0.5,
    sigma_state_prior_code = 1L,
    use_deattenuation = 1L,
    vhat_state = c(0.001, 0.002, 0.003)
  )
}

test_that("column role specifications are scalar, present, and distinct", {
  data <- make_valid_validation_data()

  expect_silent(validate_input_spec(
    data, "y", "domain", "stratum", "psu", "weight"
  ))
  expect_error(
    validate_input_spec(data, c("y", "y"), "domain", "stratum", "psu", "weight"),
    "one non-empty character"
  )
  expect_error(
    validate_input_spec(data, "", "domain", "stratum", "psu", "weight"),
    "one non-empty character"
  )
  expect_error(
    validate_input_spec(data, "y", "domain", "stratum", "psu", "psu"),
    "distinct columns"
  )
  expect_error(
    validate_input_spec(data, "missing", "domain", "stratum", "psu", "weight"),
    "Missing required variables"
  )
})

test_that("retained analysis values fail closed on invalid numerics and labels", {
  data <- make_valid_validation_data()
  expect_silent(validate_analysis_data(
    data, "y", "domain", "stratum", "psu", "weight"
  ))

  bad <- data
  bad$y <- as.character(bad$y)
  expect_error(validate_analysis_data(
    bad, "y", "domain", "stratum", "psu", "weight"
  ), "numeric or logical binary")

  bad <- data
  bad$y[1] <- Inf
  expect_error(validate_analysis_data(
    bad, "y", "domain", "stratum", "psu", "weight"
  ), "finite and non-missing")

  bad <- data
  bad$weight[1] <- 0
  expect_error(validate_analysis_data(
    bad, "y", "domain", "stratum", "psu", "weight"
  ), "strictly positive")

  bad <- data
  bad$weight[1] <- Inf
  expect_error(validate_analysis_data(
    bad, "y", "domain", "stratum", "psu", "weight"
  ), "finite and non-missing")

  bad <- data
  bad$domain[1] <- "  "
  expect_error(validate_analysis_data(
    bad, "y", "domain", "stratum", "psu", "weight"
  ), "non-blank labels")
})

test_that("preparation scalar options have exact support", {
  expect_silent(validate_preparation_options(TRUE, NULL, 0.5))
  expect_silent(validate_preparation_options(FALSE, -1.5, 1))
  expect_error(validate_preparation_options(NA, NULL, 0.5), "one non-missing logical")
  expect_error(validate_preparation_options(c(TRUE, FALSE), NULL, 0.5),
               "one non-missing logical")
  expect_error(validate_preparation_options(TRUE, Inf, 0.5), "finite numeric")
  expect_error(validate_preparation_options(TRUE, NULL, 0), "strictly positive")
  expect_error(validate_preparation_options(TRUE, NULL, Inf), "strictly positive")
})

test_that("prepare_bhf_data applies scalar-option validation", {
  data <- make_valid_validation_data()
  call_prepare <- function(...) {
    prepare_bhf_data(
      data,
      outcome = "y",
      domain = "domain",
      strata = "stratum",
      psu = "psu",
      weights = "weight",
      ...
    )
  }

  expect_error(call_prepare(use_deattenuation = NA), "one non-missing logical")
  expect_error(call_prepare(prior_alpha_mean = Inf), "finite numeric")
  expect_error(call_prepare(prior_alpha_sd = 0), "strictly positive")
})

test_that("calc_eff_n rejects empty, coercive, and invalid weights", {
  expect_equal(calc_eff_n(rep(1, 10)), 10)
  expect_lt(calc_eff_n(c(1, 1, 1, 1, 10)), 5)

  invalid <- list(
    numeric(),
    "1",
    c(1, NA_real_),
    c(1, NaN),
    c(1, Inf),
    c(1, 0),
    c(1, -1)
  )
  for (weights in invalid) {
    expect_error(calc_eff_n(weights))
  }
})

test_that("Stan-data validation has explicit dimension and support errors", {
  valid <- make_valid_stan_data()
  expect_silent(validate_stan_data(valid))

  corrupt <- function(field, value) {
    result <- valid
    result[[field]] <- value
    result
  }

  expect_error(validate_stan_data(corrupt("N", 6.5)), "positive integer")
  expect_error(validate_stan_data(corrupt("data_schema_version", 2L)),
               "integer 1")
  expect_error(validate_stan_data(corrupt("sigma_state_prior_code", 5L)),
               "integer in \\[1, 4\\]")
  expect_error(validate_stan_data(corrupt("y", valid$y[-1])), "Length of y")
  expect_error(validate_stan_data(corrupt("state_id", c(1.5, valid$state_id[-1]))),
               "integers in range")
  expect_error(validate_stan_data(corrupt("state_id", c(4L, valid$state_id[-1]))),
               "integers in range")
  expect_error(validate_stan_data(corrupt("y", c(2L, valid$y[-1]))), "binary")
  expect_error(validate_stan_data(corrupt("w_lik", c(0, valid$w_lik[-1]))),
               "strictly positive")
  expect_error(validate_stan_data(corrupt("w_lik", rep(2, valid$N))),
               "global sum-N normalization")
  expect_error(validate_stan_data(corrupt(
    "w_state_pop_share", c(0, 0.5, 0.5)
  )), "strictly positive")
  expect_error(validate_stan_data(corrupt(
    "w_state_pop_share", c(0.2, 0.3, 0.4)
  )), "sum to 1")
  expect_error(validate_stan_data(corrupt("vhat_state", c(-1, 0, 0))),
               "nonnegative")
  expect_error(validate_stan_data(corrupt("J_h", c(1L, 2L))), "summing to J")
  expect_error(validate_stan_data(corrupt("psu_start", c(1L, 2L))),
               "consecutive stratum PSU blocks")
  expect_error(validate_stan_data(corrupt(
    "psu_in_stratum_id", c(3L, valid$psu_in_stratum_id[-1])
  )), "valid within their strata")
  expect_error(validate_stan_data(corrupt(
    "psu_flat_id", c(2L, valid$psu_flat_id[-1])
  )), "stratum-block reconstruction")
})
