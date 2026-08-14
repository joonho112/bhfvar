test_that("Phase 7 scientific tests pin the frozen oracle and tolerance", {
  paths <- c(
    "estimand_oracle.R" = testthat::test_path(
      "..", "oracle", "estimand_oracle.R"
    ),
    "smoke_oracle.R" = testthat::test_path(
      "..", "oracle", "smoke_oracle.R"
    ),
    "manifest.txt" = testthat::test_path("..", "oracle", "manifest.txt"),
    "test-estimand-oracle.R" = testthat::test_path(
      "test-estimand-oracle.R"
    )
  )
  observed <- vapply(paths, .bhf_sha256_file, character(1))
  expect_identical(observed, .bhf_science_frozen_hashes)

  manifest <- readLines(paths[["manifest.txt"]], warn = FALSE)
  expect_true("oracle_version: 1.0.0" %in% manifest)
  expect_true(paste0(
    "default_tolerance: abs(x-y)<=1e-10+1e-8*max(abs(x),abs(y))"
  ) %in% manifest)
  expect_identical(
    .bhf_science_tolerance,
    c(absolute = 1e-10, relative = 1e-8)
  )
})

test_that("seeded randomized A, A-star, B, latent, and gaps close algebraically", {
  fixture <- .bhf_science_random_fixture()
  result <- .bhf_science_oracle(
    fixture$design,
    fixture$draws,
    vhat = fixture$vhat,
    return_individual = TRUE
  )

  numeric_surfaces <- list(
    latent = result$latent,
    A_state = result$A$p_state,
    A = result$A$summary,
    A_star = result$A_star$summary,
    B_state = result$B$p_state,
    B_individual = result$B$p_individual,
    B_binomial_state = result$B$within_bernoulli_state,
    B_mixture_state = result$B$within_mixture_state,
    B = result$B$summary,
    B_minus_A = result$gaps$B_minus_A,
    A_minus_A_star = result$gaps$A_minus_A_star
  )
  expect_true(all(vapply(
    numeric_surfaces,
    function(x) all(is.finite(x)),
    logical(1)
  )))
  expect_true(all(result$A$p_state >= 0 & result$A$p_state <= 1))
  expect_true(all(result$B$p_state >= 0 & result$B$p_state <= 1))
  expect_true(all(result$B$p_individual >= 0 & result$B$p_individual <= 1))

  for (summary in list(
    A = result$A$summary,
    A_star = result$A_star$summary,
    B = result$B$summary
  )) {
    expect_true(all(summary[, c("between", "within", "total")] >= 0))
    expect_true(all(summary[, "proportion"] >= 0 &
                      summary[, "proportion"] <= 1))
    .bhf_science_expect_close(
      summary[, "total"],
      summary[, "between"] + summary[, "within"],
      "total equals between plus within"
    )
  }
  for (summary in list(A = result$A$summary, B = result$B$summary)) {
    .bhf_science_expect_close(
      summary[, "total"],
      summary[, "mean"] * (1 - summary[, "mean"]),
      "binary total equals p-bar times one-minus-p-bar"
    )
  }

  .bhf_science_expect_close(
    result$B$within_bernoulli_state + result$B$within_mixture_state,
    result$B$p_state * (1 - result$B$p_state),
    "B state binomial plus mixture identity"
  )
  .bhf_science_expect_close(
    result$A_star$summary[, "mean"], result$A$summary[, "mean"],
    "A-star preserves A mean"
  )
  .bhf_science_expect_close(
    result$A_star$summary[, "within"], result$A$summary[, "within"],
    "A-star preserves A within"
  )
  expect_true(all(result$A_star$summary[, "between"] <=
                    result$A$summary[, "between"] + 1e-15))

  .bhf_science_expect_close(
    result$gaps$B_minus_A,
    result$B$summary - result$A$summary,
    "signed B-minus-A gap"
  )
  .bhf_science_expect_close(
    result$gaps$A_minus_A_star,
    result$A$summary - result$A_star$summary,
    "signed A-minus-A-star gap"
  )

  latent <- result$latent
  .bhf_science_expect_close(
    latent[, "total"],
    latent[, "var_state"] + latent[, "var_stratum"] +
      latent[, "var_psu"] + latent[, "var_logistic"],
    "latent total identity"
  )
  expect_true(all(latent[, c(
    "var_state", "var_stratum", "var_psu", "var_logistic", "total"
  )] >= 0))
  expect_true(all(latent[, c(
    "icc_state", "icc_stratum", "icc_psu", "icc_logistic"
  )] >= 0))
  .bhf_science_expect_close(
    rowSums(latent[, c(
      "icc_state", "icc_stratum", "icc_psu", "icc_logistic"
    )]),
    rep(1, nrow(latent)),
    "latent ICCs including level one sum to one"
  )
})

test_that("the randomized scientific fixture is reproducible and RNG-local", {
  set.seed(909L)
  expected_next <- runif(1)
  set.seed(909L)
  first <- .bhf_science_random_fixture(.bhf_science_seed)
  observed_next <- runif(1)
  second <- .bhf_science_random_fixture(.bhf_science_seed)

  expect_identical(observed_next, expected_next)
  expect_identical(first, second)
  expect_identical(first$draws$draw_id, paste0("d", 1:32))
  expect_equal(sum(first$design$w_lik), length(first$design$w_lik),
               tolerance = 1e-12)
})
