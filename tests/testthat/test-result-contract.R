test_that("shared interval contract has dynamic stable labels", {
  interval_95 <- .bhf_interval_contract(0.95)
  interval_90 <- .bhf_interval_contract(0.90)

  expect_identical(unname(interval_95$labels), c("2.5%", "50%", "97.5%"))
  expect_identical(unname(interval_90$labels), c("5%", "50%", "95%"))
  expect_identical(interval_95$interval_label, "95% credible interval")
  expect_identical(interval_90$interval_label, "90% credible interval")
  expect_named(
    .bhf_draw_summary(c(0, 0, 0), interval_90),
    c("mean", "sd", "lower", "median", "upper", "prob")
  )
  expect_true(all(is.finite(.bhf_draw_summary(c(0, 0), interval_95))))

  for (value in list(0, 1, -0.1, 1.1, NA_real_, Inf, c(0.9, 0.95), "0.9")) {
    expect_error(.bhf_interval_contract(value), class = "bhf_interval_error")
  }
})

test_that("variance API separates latent, A, A-star, B, and signed gaps", {
  local_result_draw_backend()
  fit <- make_result_contract_fit(TRUE)
  result <- variance_decomposition(fit, prob = 0.90)

  expect_s3_class(result, "bhf_variance_decomposition")
  expect_named(
    result,
    c("schema_version", "interval", "latent", "A", "A_star", "B",
      "gaps", "summary_table")
  )
  expect_identical(result$schema_version, "0.5.0")
  expect_identical(result$interval$labels,
                   c(lower = "5%", median = "50%", upper = "95%"))
  expect_false(any(c("logit", "prob", "deatten") %in% names(result)))
  expect_true(result$A_star$available)
  expect_identical(result$A_star$provenance$mode, "supplied")
  expect_named(
    result$summary_table,
    c("section", "estimand", "scale", "component", "mean", "sd",
      "lower", "median", "upper", "prob", "available")
  )
  expect_true(all(result$summary_table$prob == 0.90))
  expect_true(all(result$summary_table$available))

  draws <- fit$.draws
  design_mean <- result$gaps$B_minus_A$mean[
    result$gaps$B_minus_A$component == "between"
  ]
  noise_mean <- result$gaps$A_minus_A_star$mean[
    result$gaps$A_minus_A_star$component == "between"
  ]
  expect_equal(design_mean, mean(draws$var_between_B - draws$var_between_A))
  expect_equal(
    noise_mean,
    mean(draws$var_between_A - draws$var_between_A_star)
  )
})

test_that("disabled A-star is unavailable without reinterpreting zero carrier", {
  local_result_draw_backend()
  fit <- make_result_contract_fit(FALSE)
  fit$.draws$var_total_A <- rep(0, 5)
  fit$.draws$prop_between_A <- rep(0, 5)
  fit$.draws$A_proportion_defined <- rep(0L, 5)
  fit$.draws$var_total_B <- rep(0, 5)
  fit$.draws$prop_between_B <- rep(0, 5)
  fit$.draws$B_proportion_defined <- rep(0L, 5)

  result <- variance_decomposition(fit)
  expect_false(result$A_star$available)
  expect_null(result$A_star$summary)
  expect_null(result$A_star$correction)
  expect_null(result$A_star$flags)
  expect_null(result$gaps$A_minus_A_star)
  unavailable <- subset(result$summary_table, estimand == "A_star")
  expect_true(all(!unavailable$available))
  expect_true(all(is.na(unavailable$mean)))
  available <- subset(result$summary_table, available)
  expect_false(any(is.nan(available$mean)))
  expect_equal(
    subset(result$A$summary, component == "proportion")$mean,
    0
  )
  expect_false(result$A$proportion_defined$any)
  expect_false(result$B$proportion_defined$any)
})

test_that("domain estimates use explicit A/B vocabulary and canonical order", {
  local_result_draw_backend()
  fit <- make_result_contract_fit(TRUE)

  a <- domain_estimates(fit, estimand = "A", prob = 0.90)
  b <- domain_estimates(fit, estimand = "B", prob = 0.95)
  expected_names <- c(
    "domain", "domain_id", "estimand", "mean", "sd", "lower", "median",
    "upper", "prob", "pop_share", "n"
  )
  expect_named(a, expected_names)
  expect_named(b, expected_names)
  expect_identical(a$domain, c("A", "B", "C"))
  expect_identical(a$domain_id, 1:3)
  expect_identical(a$estimand, rep("A", 3))
  expect_identical(b$estimand, rep("B", 3))
  expect_equal(a$mean, colMeans(fit$.draws$p_state_A))
  expect_equal(b$mean, colMeans(fit$.draws$p_state_B))
  expect_identical(a$n, c(2L, 2L, 2L))
  expect_equal(a$pop_share, c(0.5, 0.3, 0.2))
  expect_identical(attr(a, "schema_version"), "0.5.0")
  expect_identical(attr(a, "interval")$labels,
                   c(lower = "5%", median = "50%", upper = "95%"))

  expect_warning(
    legacy_a <- domain_estimates(fit, type = "conditional"),
    class = "bhf_deprecated_argument_warning"
  )
  expect_identical(legacy_a$estimand, rep("A", 3))
  expect_error(domain_estimates(fit, type = "marginal"),
               class = "bhf_legacy_marginal_error")
  expect_error(domain_estimates(fit, estimand = "A", type = "conditional"),
               class = "bhf_argument_error")
})

test_that("overall estimates verify population-share reconstruction", {
  local_result_draw_backend()
  fit <- make_result_contract_fit(TRUE)

  a <- overall_estimate(fit, "A", prob = 0.90)
  b <- overall_estimate(fit, "B", prob = 0.95)
  expect_s3_class(a, "bhf_overall_estimate")
  expect_named(
    a,
    c("schema_version", "estimand", "scale", "interval",
      "population_shares", "mean", "sd", "lower", "median", "upper",
      "prob")
  )
  expect_equal(a$mean, mean(fit$.draws$p_bar_A))
  expect_equal(b$mean, mean(fit$.draws$p_bar_B))
  expect_identical(a$estimand, "A")
  expect_identical(b$estimand, "B")
  expect_equal(a$population_shares$values, c(A = 0.5, B = 0.3, C = 0.2))
  expect_identical(
    a$population_shares$provenance$source, "external_known"
  )

  broken <- fit
  broken$.draws$p_bar_B <- broken$.draws$p_bar_B + 0.01
  expect_error(overall_estimate(broken, "B"),
               class = "bhf_result_invariant_error")
})

test_that("log-likelihood API distinguishes kind and aggregation scope", {
  local_result_draw_backend()
  fit <- make_result_contract_fit(TRUE)
  raw <- fit$.draws$log_lik_raw
  pseudo <- fit$.draws$log_lik_pseudo

  expect_warning(
    observation <- log_lik(fit, kind = "pseudo", aggregate = "observation"),
    class = "bhf_pseudo_loo_warning"
  )
  expect_warning(
    raw_observation <- log_lik(fit, kind = "raw"),
    class = "bhf_raw_log_lik_warning"
  )
  expect_s3_class(observation, "bhf_log_lik")
  expect_equal(unname(observation), unname(pseudo), ignore_attr = TRUE)
  expect_equal(unname(raw_observation), unname(raw), ignore_attr = TRUE)
  expect_identical(colnames(observation), c("2", "4", "5", "8", "9", "10"))
  expect_identical(attr(observation, "kind"), "pseudo")
  expect_false(attr(observation, "scope")$ordinary_loo_supported)
  expect_false(attr(observation, "scope")$cluster_loo_supported)

  expect_warning(psu <- log_lik(fit, aggregate = "psu"),
                 class = "bhf_pseudo_loo_warning")
  expect_warning(stratum <- log_lik(fit, aggregate = "stratum"),
                 class = "bhf_pseudo_loo_warning")
  expected_psu <- cbind(
    rowSums(pseudo[, 1:2, drop = FALSE]), pseudo[, 3], pseudo[, 4],
    rowSums(pseudo[, 5:6, drop = FALSE])
  )
  expected_stratum <- cbind(
    rowSums(pseudo[, 1:3, drop = FALSE]),
    rowSums(pseudo[, 4:6, drop = FALSE])
  )
  expect_equal(unname(psu), unname(expected_psu), ignore_attr = TRUE)
  expect_equal(unname(stratum), unname(expected_stratum), ignore_attr = TRUE)
  expect_identical(colnames(psu), c("H1::P1", "H1::P2", "H2::P1", "H2::P2"))
  expect_identical(colnames(stratum), c("H1", "H2"))
  expect_identical(attr(psu, "aggregation"), "psu")
  expect_match(attr(psu, "scope")$caveat, "not leave-one-unit refitting")
})

test_that("posterior extraction boundary fails closed on missing fields", {
  local_result_draw_backend()
  fit <- make_result_contract_fit(TRUE)
  fit$.draws$p_state_A <- NULL
  expect_error(domain_estimates(fit, "A"), class = "bhf_draw_contract_error")

  fit <- make_result_contract_fit(TRUE)
  fit$.draws$var_between_A <- matrix(1:10, ncol = 2)
  expect_error(variance_decomposition(fit), class = "bhf_draw_shape_error")
})
