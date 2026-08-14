test_that("result extractors require a versioned bhf_fit object", {
  functions <- list(
    variance_decomposition = variance_decomposition,
    domain_estimates = domain_estimates,
    overall_estimate = overall_estimate,
    log_lik = log_lik
  )
  for (fun in functions) {
    expect_error(fun(list(x = 1)), class = "bhf_fit_class_error")
  }

  legacy <- structure(
    list(stanfit = list(), data = list(stan_data = list())),
    class = c("bhf_fit", "list")
  )
  for (fun in functions) {
    expect_error(fun(legacy), class = "bhf_unsupported_schema_error")
  }
})

test_that("domain and log-likelihood choices validate before fit class", {
  expect_error(
    domain_estimates(list(), estimand = "invalid"),
    class = "bhf_argument_error"
  )
  expect_error(
    domain_estimates(list(), type = "invalid"),
    class = "bhf_deprecated_argument_error"
  )
  expect_error(
    domain_estimates(list(), type = "marginal"),
    class = "bhf_legacy_marginal_error"
  )
  expect_error(log_lik(list(), kind = "invalid"),
               class = "bhf_argument_error")
  expect_error(log_lik(list(), aggregate = "cluster"),
               class = "bhf_argument_error")
})
