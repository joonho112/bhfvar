test_that("weighted state, stratum, and every PSU block center in every draw", {
  fixture <- .bhf_science_random_fixture(.bhf_science_seed + 1L, draws = 48L)
  residual <- .bhf_science_assert_centered(fixture)

  expect_identical(length(residual$state), 48L)
  expect_identical(length(residual$stratum), 48L)
  expect_identical(dim(residual$psu), c(48L, 3L))
  expect_lt(max(abs(residual$state)), 1e-12)
  expect_lt(max(abs(residual$stratum)), 1e-12)
  expect_lt(max(abs(residual$psu)), 1e-12)
  expect_identical(unname(residual$psu[, 3L]), rep(0, 48L))
})

test_that("centering regression failures identify component, draw, and block", {
  fixture <- .bhf_science_random_fixture(.bhf_science_seed + 2L, draws = 6L)

  bad_state <- fixture
  bad_state$draws$u_state[3L, 1L] <-
    bad_state$draws$u_state[3L, 1L] + 0.1
  expect_error(
    .bhf_science_assert_centered(bad_state),
    "state centering failed at draw d3"
  )

  bad_stratum <- fixture
  bad_stratum$draws$u_stratum[4L, 2L] <-
    bad_stratum$draws$u_stratum[4L, 2L] - 0.2
  expect_error(
    .bhf_science_assert_centered(bad_stratum),
    "stratum centering failed at draw d4"
  )

  bad_psu <- fixture
  second_block_first_psu <- fixture$psu_start[2L]
  bad_psu$draws$u_psu[5L, second_block_first_psu] <-
    bad_psu$draws$u_psu[5L, second_block_first_psu] + 0.3
  expect_error(
    .bhf_science_assert_centered(bad_psu),
    "PSU centering failed at draw d5, stratum H2"
  )
})
