oracle_path <- testthat::test_path("..", "oracle", "estimand_oracle.R")
sys.source(oracle_path, envir = environment())

make_exact_estimand_oracle_fixture <- function(informative = FALSE) {
  L <- log(3)
  design <- list(
    row_id = sprintf("r%d", 1:8),
    domain_labels = c("A", "B"),
    stratum_labels = c("H1", "H2"),
    psu_labels = c("J1", "J2", "J3", "J4"),
    state_id = c(rep(1L, 4), rep(2L, 4)),
    stratum_id = rep(c(1L, 1L, 2L, 2L), 2),
    psu_flat_id = rep(1:4, 2),
    w_lik = if (informative) {
      c(2.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 2.5)
    } else rep(1, 8),
    population_share = c(A = 0.5, B = 0.5)
  )
  draws <- list(
    draw_id = "d1",
    alpha = 0,
    u_state = matrix(c(-L, L), 1, dimnames = list("d1", c("A", "B"))),
    u_stratum = matrix(
      if (informative) c(L / 2, -L / 2) else c(0, 0),
      1, dimnames = list("d1", c("H1", "H2"))
    ),
    u_psu = matrix(
      if (informative) c(L / 2, -L / 2, L / 2, -L / 2) else rep(0, 4),
      1, dimnames = list("d1", c("J1", "J2", "J3", "J4"))
    ),
    sigma_state = L,
    sigma_stratum = if (informative) L / 2 else 0,
    sigma_psu = if (informative) L / 2 else 0
  )
  list(design = design, draws = draws)
}

test_that("no-design gold fixture gives exact A=A*=B", {
  fixture <- make_exact_estimand_oracle_fixture(FALSE)
  result <- bhf_reference_oracle(
    fixture$design, fixture$draws, vhat = c(A = 0, B = 0),
    return_individual = TRUE
  )
  expected <- c(
    mean = 1 / 2, between = 1 / 16, within = 3 / 16,
    total = 1 / 4, proportion = 1 / 4
  )

  expect_equal(unname(result$A$p_state[1, ]), c(1 / 4, 3 / 4),
               tolerance = 1e-12)
  expect_equal(result$A$summary[1, ], expected, tolerance = 1e-12)
  expect_equal(result$B$p_state, result$A$p_state, tolerance = 1e-12)
  expect_equal(result$B$summary, result$A$summary, tolerance = 1e-12)
  expect_equal(unname(result$B$within_mixture), 0, tolerance = 1e-12)
  expect_equal(result$A_star$summary, result$A$summary, tolerance = 1e-12)
  expect_equal(result$gaps$B_minus_A, matrix(0, 1, 5), tolerance = 1e-12,
               ignore_attr = TRUE)

  L <- log(3)
  expected_total <- L^2 + pi^2 / 3
  expect_equal(result$latent[1, "total"], expected_total, tolerance = 1e-12)
  expect_equal(result$latent[1, "icc_state"], L^2 / expected_total,
               tolerance = 1e-12)
})

test_that("informative-design gold fixture matches exact B components", {
  fixture <- make_exact_estimand_oracle_fixture(TRUE)
  result <- bhf_reference_oracle(fixture$design, fixture$draws)

  expect_equal(unname(result$B$p_state[1, ]), c(31 / 80, 49 / 80),
               tolerance = 1e-12)
  expect_equal(unname(result$B$within_bernoulli), 343 / 1600,
               tolerance = 1e-12)
  expect_equal(unname(result$B$within_mixture), 147 / 6400,
               tolerance = 1e-12)
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
  expect_false(result$A_star$available)
  expect_null(result$A_star$summary)
  expect_null(result$gaps$A_minus_A_star)
  expect_equal(
    result$B$within_bernoulli_state + result$B$within_mixture_state,
    result$B$p_state * (1 - result$B$p_state),
    tolerance = 1e-12
  )
  expect_equal(
    result$B$summary[1, "total"],
    result$B$summary[1, "mean"] * (1 - result$B$summary[1, "mean"]),
    tolerance = 1e-12
  )
})

test_that("A-star zero, boundary, and overshoot rules are exact", {
  fixture <- make_exact_estimand_oracle_fixture(FALSE)
  zero <- bhf_reference_oracle(
    fixture$design, fixture$draws, c(A = 0, B = 0)
  )
  boundary <- bhf_reference_oracle(
    fixture$design, fixture$draws, c(B = 3 / 32, A = 1 / 32)
  )
  overshoot <- bhf_reference_oracle(
    fixture$design, fixture$draws, c(A = 1 / 16, B = 3 / 16)
  )

  expect_equal(zero$A_star$summary, zero$A$summary, tolerance = 1e-12)
  expect_equal(boundary$A_star$summary[1, ], c(
    mean = 1 / 2, between = 0, within = 3 / 16,
    total = 3 / 16, proportion = 0
  ), tolerance = 1e-12)
  expect_equal(overshoot$A_star$summary[1, "between"], 0,
               tolerance = 1e-12)
  expect_true(boundary$A_star$at_boundary[1])
  expect_false(boundary$A_star$truncated[1])
  expect_true(overshoot$A_star$truncated[1])
  expect_equal(boundary$gaps$A_minus_A_star[1, ], c(
    mean = 0, between = 1 / 16, within = 0,
    total = 1 / 16, proportion = 1 / 4
  ), tolerance = 1e-12)
})

test_that("multi-draw output stacks single-draw results without dropping", {
  fixture <- make_exact_estimand_oracle_fixture(TRUE)
  two <- fixture$draws
  two$draw_id <- c("d1", "d2")
  two$alpha <- rep(two$alpha, 2)
  two$u_state <- rbind(two$u_state, two$u_state)
  rownames(two$u_state) <- two$draw_id
  two$u_stratum <- rbind(two$u_stratum, two$u_stratum)
  rownames(two$u_stratum) <- two$draw_id
  two$u_psu <- rbind(two$u_psu, two$u_psu)
  rownames(two$u_psu) <- two$draw_id
  two$sigma_state <- rep(two$sigma_state, 2)
  two$sigma_stratum <- rep(two$sigma_stratum, 2)
  two$sigma_psu <- rep(two$sigma_psu, 2)

  stacked <- bhf_reference_oracle(fixture$design, two)
  single <- bhf_reference_oracle(fixture$design, fixture$draws)
  expect_identical(dim(stacked$A$p_state), c(2L, 2L))
  expect_equal(stacked$A$summary[1, ], single$A$summary[1, ],
               tolerance = 1e-12)
  expect_equal(stacked$A$summary[2, ], single$A$summary[1, ],
               tolerance = 1e-12)
  expect_equal(stacked$B$summary[1, ], stacked$B$summary[2, ],
               tolerance = 1e-12)
})

test_that("row and exact-name permutations preserve state summaries", {
  fixture <- make_exact_estimand_oracle_fixture(TRUE)
  reference <- bhf_reference_oracle(fixture$design, fixture$draws,
                                    c(A = 1 / 32, B = 3 / 32))

  permutation <- c(8L, 1L, 6L, 3L, 5L, 2L, 7L, 4L)
  row_fields <- c("row_id", "state_id", "stratum_id", "psu_flat_id", "w_lik")
  for (field in row_fields) {
    fixture$design[[field]] <- fixture$design[[field]][permutation]
  }
  fixture$design$population_share <- c(B = 0.5, A = 0.5)
  fixture$draws$u_state <- fixture$draws$u_state[, c("B", "A"), drop = FALSE]

  observed <- bhf_reference_oracle(
    fixture$design, fixture$draws, c(B = 3 / 32, A = 1 / 32)
  )
  expect_equal(observed$A$summary, reference$A$summary, tolerance = 1e-12)
  expect_equal(observed$A_star$summary, reference$A_star$summary,
               tolerance = 1e-12)
  expect_equal(observed$B$summary, reference$B$summary, tolerance = 1e-12)
})

test_that("total-zero and invalid-universe policies fail or mark undefined", {
  fixture <- make_exact_estimand_oracle_fixture(FALSE)
  fixture$draws$alpha <- -1000
  fixture$draws$u_state[,] <- 0
  zero <- bhf_reference_oracle(fixture$design, fixture$draws)
  zero_star <- bhf_reference_oracle(
    fixture$design, fixture$draws, vhat = c(A = 0, B = 0)
  )
  expect_equal(zero$A$summary[1, "total"], 0)
  expect_equal(zero$A$summary[1, "proportion"], 0)
  expect_false(zero$A$proportion_defined[1])
  expect_equal(zero$B$summary[1, "total"], 0)
  expect_equal(zero$B$summary[1, "proportion"], 0)
  expect_false(zero$B$proportion_defined[1])
  expect_equal(zero_star$A_star$summary[1, "total"], 0)
  expect_equal(zero_star$A_star$summary[1, "proportion"], 0)
  expect_false(zero_star$A_star$proportion_defined[1])

  bad <- make_exact_estimand_oracle_fixture(FALSE)
  bad$design$domain_labels <- c("A", "B", "C")
  bad$design$population_share <- c(A = 0.4, B = 0.4, C = 0.2)
  expect_error(bhf_reference_oracle(bad$design, bad$draws),
               "complete consecutive observed levels")
})

test_that("oracle fails closed on dimensions, labels, and weight support", {
  fixture <- make_exact_estimand_oracle_fixture(FALSE)

  bad <- fixture
  colnames(bad$draws$u_state) <- NULL
  expect_error(bhf_reference_oracle(bad$design, bad$draws),
               "draw-by-level matrix")

  bad <- fixture
  rownames(bad$draws$u_state) <- "wrong"
  expect_error(bhf_reference_oracle(bad$design, bad$draws),
               "draw-by-level matrix")

  bad <- fixture
  bad$design$psu_flat_id[c(5, 7)] <- c(3L, 1L)
  expect_error(bhf_reference_oracle(bad$design, bad$draws),
               "nested in exactly one stratum")

  bad <- fixture
  bad$design$w_lik[1] <- 2
  expect_error(bhf_reference_oracle(bad$design, bad$draws),
               "globally sum to N")

  bad <- fixture
  bad$design$population_share <- c(A = 0.5, C = 0.5)
  expect_error(bhf_reference_oracle(bad$design, bad$draws),
               "exact-set named vector")

  bad <- fixture
  bad$draws$sigma_state <- -1
  expect_error(bhf_reference_oracle(bad$design, bad$draws),
               "finite and nonnegative")
})

test_that("oracle source has no package dependency calls", {
  source_text <- paste(readLines(oracle_path, warn = FALSE), collapse = "\n")
  forbidden <- c(
    "library\\s*\\(", "require\\s*\\(", "requireNamespace\\s*\\(",
    "loadNamespace\\s*\\(", "getFromNamespace\\s*\\(", "::", ":::"
  )
  for (pattern in forbidden) {
    expect_false(grepl(pattern, source_text, perl = TRUE))
  }
  fixture <- make_exact_estimand_oracle_fixture(FALSE)
  result <- bhf_reference_oracle(fixture$design, fixture$draws)
  expect_identical(result$schema_version, "1.0.0")
  expect_identical(
    result$formula_authority,
    "math14030512-probability-decomposition-v1"
  )
  expect_identical(result$dimensions,
                   c(draws = 1L, observations = 8L, states = 2L,
                     strata = 2L, psus = 4L))
})
