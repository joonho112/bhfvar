test_that("disabled, zero, equal, and overshoot vhat follow the frozen A-star policy", {
  fixture <- .bhf_science_exact_fixture(FALSE)
  disabled <- .bhf_science_oracle(fixture$design, fixture$draws)
  zero <- .bhf_science_oracle(
    fixture$design, fixture$draws, vhat = c(A = 0, B = 0)
  )
  equal <- .bhf_science_oracle(
    fixture$design, fixture$draws,
    vhat = c(B = 3 / 32, A = 1 / 32)
  )
  overshoot <- .bhf_science_oracle(
    fixture$design, fixture$draws,
    vhat = c(A = 1 / 16, B = 3 / 16)
  )

  expect_false(disabled$A_star$available)
  expect_null(disabled$A_star$summary)
  expect_null(disabled$A_star$correction)
  expect_null(disabled$gaps$A_minus_A_star)

  expect_true(zero$A_star$available)
  expect_equal(zero$A_star$summary, zero$A$summary, tolerance = 1e-12)
  expect_equal(zero$A_star$correction, 0, tolerance = 0)
  expect_false(zero$A_star$at_boundary[1])
  expect_false(zero$A_star$truncated[1])

  expect_equal(equal$A_star$correction, 1 / 16, tolerance = 1e-12)
  expect_equal(equal$A_star$summary[1, "between"], 0, tolerance = 0)
  expect_equal(equal$A_star$summary[1, "within"], 3 / 16,
               tolerance = 1e-12)
  expect_equal(equal$A_star$summary[1, "total"], 3 / 16,
               tolerance = 1e-12)
  expect_equal(equal$A_star$summary[1, "proportion"], 0, tolerance = 0)
  expect_true(equal$A_star$at_boundary[1])
  expect_false(equal$A_star$truncated[1])

  expect_equal(overshoot$A_star$correction, 1 / 8, tolerance = 1e-12)
  expect_equal(overshoot$A_star$summary[1, "between"], 0, tolerance = 0)
  expect_true(overshoot$A_star$at_boundary[1])
  expect_true(overshoot$A_star$truncated[1])
  expect_true(all(overshoot$A_star$summary >= 0))
})

test_that("named vhat permutations preserve A-star and invalid inputs fail closed", {
  fixture <- .bhf_science_exact_fixture(FALSE)
  canonical <- .bhf_science_oracle(
    fixture$design, fixture$draws,
    vhat = c(A = 1 / 32, B = 3 / 32)
  )
  permuted <- .bhf_science_oracle(
    fixture$design, fixture$draws,
    vhat = c(B = 3 / 32, A = 1 / 32)
  )
  expect_identical(permuted$A_star, canonical$A_star)
  expect_identical(permuted$gaps$A_minus_A_star,
                   canonical$gaps$A_minus_A_star)

  invalid <- list(
    c(A = -0.01, B = 0.02),
    c(0.01, 0.02),
    c(A = 0.01),
    c(A = 0.01, C = 0.02),
    stats::setNames(c(0.01, 0.02), c("A", "A")),
    c(A = Inf, B = 0.02),
    c(A = NA_real_, B = 0.02)
  )
  for (vhat in invalid) {
    expect_error(
      .bhf_science_oracle(fixture$design, fixture$draws, vhat = vhat),
      "vhat must be a finite exact-set named vector"
    )
  }
})

test_that("package supplied-vhat resolver shares the boundary input contract", {
  labels <- c("A", "B")
  valid <- resolve_sampling_variances(
    deattenuation = "supplied",
    sampling_variances = c(B = 3 / 32, A = 1 / 32),
    sampling_variance_method = "external_taylor",
    domain_labels = labels
  )
  expect_identical(valid$named_values, c(A = 1 / 32, B = 3 / 32))
  expect_identical(valid$stan_values, c(1 / 32, 3 / 32))

  for (vhat in list(
    c(A = -0.01, B = 0.02),
    c(0.01, 0.02),
    c(A = 0.01),
    c(A = 0.01, C = 0.02),
    c(A = Inf, B = 0.02)
  )) {
    expect_error(resolve_sampling_variances(
      deattenuation = "supplied",
      sampling_variances = vhat,
      sampling_variance_method = "external_taylor",
      domain_labels = labels
    ))
  }
})
