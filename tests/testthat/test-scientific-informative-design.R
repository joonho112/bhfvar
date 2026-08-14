test_that("informative design has the frozen B-minus-A direction and magnitude", {
  fixture <- .bhf_science_exact_fixture(TRUE)
  result <- .bhf_science_oracle(fixture$design, fixture$draws)

  expect_equal(
    unname(result$A$p_state[1, ]),
    c(1 / 4, 3 / 4),
    tolerance = 1e-12
  )
  expect_equal(
    unname(result$B$p_state[1, ]),
    c(31 / 80, 49 / 80),
    tolerance = 1e-12
  )
  expect_equal(
    result$A$summary[1, ],
    c(
      mean = 1 / 2, between = 1 / 16, within = 3 / 16,
      total = 1 / 4, proportion = 1 / 4
    ),
    tolerance = 1e-12
  )
  expect_equal(
    result$B$summary[1, ],
    c(
      mean = 1 / 2, between = 81 / 6400, within = 1519 / 6400,
      total = 1 / 4, proportion = 81 / 1600
    ),
    tolerance = 1e-12
  )
  expect_equal(
    result$gaps$B_minus_A[1, ],
    c(
      mean = 0, between = -319 / 6400, within = 319 / 6400,
      total = 0, proportion = -319 / 1600
    ),
    tolerance = 1e-12
  )
  expect_equal(unname(result$B$within_bernoulli), 343 / 1600,
               tolerance = 1e-12)
  expect_equal(unname(result$B$within_mixture), 147 / 6400,
               tolerance = 1e-12)

  gap <- result$gaps$B_minus_A[1, ]
  expect_lt(gap[["between"]], 0)
  expect_gt(gap[["within"]], 0)
  expect_lt(gap[["proportion"]], 0)
  expect_equal(gap[["mean"]], 0, tolerance = 1e-12)
  expect_equal(gap[["total"]], 0, tolerance = 1e-12)
  expect_gt(result$B$p_state[1, "A"], result$A$p_state[1, "A"])
  expect_lt(result$B$p_state[1, "B"], result$A$p_state[1, "B"])
})

test_that("informative-design labels encode the approved signed interpretation", {
  fixture <- .bhf_science_exact_fixture(TRUE)
  result <- .bhf_science_oracle(fixture$design, fixture$draws)

  expect_identical(names(result$gaps), c("B_minus_A", "A_minus_A_star"))
  expect_identical(
    colnames(result$gaps$B_minus_A),
    c("mean", "between", "within", "total", "proportion")
  )
  expect_null(result$gaps$A_minus_A_star)
  expect_false(result$A_star$available)
})
