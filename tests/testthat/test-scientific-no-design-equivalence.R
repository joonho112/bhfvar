test_that("zero design effects make independent A and B formulas equivalent", {
  fixture <- .bhf_science_random_fixture(.bhf_science_seed + 3L, draws = 24L)
  fixture$draws$u_stratum[,] <- 0
  fixture$draws$u_psu[,] <- 0
  fixture$draws$sigma_stratum[] <- 0
  fixture$draws$sigma_psu[] <- 0
  zero_vhat <- stats::setNames(
    numeric(length(fixture$design$domain_labels)),
    fixture$design$domain_labels
  )
  result <- .bhf_science_oracle(
    fixture$design,
    fixture$draws,
    vhat = zero_vhat,
    return_individual = TRUE
  )

  direct_A <- t(vapply(seq_along(fixture$draws$alpha), function(draw) {
    stats::plogis(fixture$draws$alpha[draw] +
                    fixture$draws$u_state[draw, ])
  }, numeric(length(fixture$design$domain_labels))))
  dimnames(direct_A) <- dimnames(result$A$p_state)

  .bhf_science_expect_close(result$A$p_state, direct_A,
                            "independent direct A probabilities")
  .bhf_science_expect_close(result$B$p_state, direct_A,
                            "zero-design B state means equal A")
  .bhf_science_expect_close(result$B$summary, result$A$summary,
                            "zero-design decomposition A equals B")
  .bhf_science_expect_close(
    result$gaps$B_minus_A,
    matrix(0, nrow(result$gaps$B_minus_A), ncol(result$gaps$B_minus_A)),
    "zero-design B-minus-A gap"
  )
  .bhf_science_expect_close(
    result$B$within_mixture,
    rep(0, length(result$B$within_mixture)),
    "zero-design within-state mixture is zero"
  )
  .bhf_science_expect_close(result$A_star$summary, result$A$summary,
                            "zero-vhat A-star equals A")
})

test_that("the exact zero-design gold case remains rational", {
  fixture <- .bhf_science_exact_fixture(FALSE)
  result <- .bhf_science_oracle(
    fixture$design,
    fixture$draws,
    vhat = c(A = 0, B = 0)
  )
  expected <- c(
    mean = 1 / 2,
    between = 1 / 16,
    within = 3 / 16,
    total = 1 / 4,
    proportion = 1 / 4
  )

  expect_equal(result$A$summary[1, ], expected, tolerance = 1e-12)
  expect_equal(result$B$summary[1, ], expected, tolerance = 1e-12)
  expect_equal(result$A_star$summary[1, ], expected, tolerance = 1e-12)
})
