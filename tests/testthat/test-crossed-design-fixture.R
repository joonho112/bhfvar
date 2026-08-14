test_that("tiny crossed-design fixture is deterministic and complete", {
  fixture <- make_tiny_crossed_design_fixture()

  expect_identical(fixture, make_tiny_crossed_design_fixture())
  expect_named(fixture, c("data", "known", "truth"))
  expect_identical(
    fixture$truth$dimensions,
    c(N = 12L, S = 3L, H = 2L, J = 4L)
  )
  expect_identical(unique(fixture$data$state), fixture$truth$state_levels)
  expect_identical(unique(fixture$data$stratum), fixture$truth$stratum_levels)
})

test_that("states are crossed with strata and PSUs are nested by composite key", {
  fixture <- make_tiny_crossed_design_fixture()
  data <- fixture$data
  keys <- unname(fixture$truth$composite_psu_key)

  state_by_stratum <- table(data$state, data$stratum)
  expect_true(all(state_by_stratum > 0))

  strata_per_raw_psu <- vapply(
    split(data$stratum, data$psu),
    function(x) length(unique(x)),
    integer(1)
  )
  expect_identical(unname(strata_per_raw_psu), c(2L, 2L))

  strata_per_composite_psu <- vapply(
    split(data$stratum, keys),
    function(x) length(unique(x)),
    integer(1)
  )
  expect_true(all(strata_per_composite_psu == 1L))
  expect_identical(unique(keys), fixture$truth$composite_psu_levels)
})

test_that("crossed-design truth contains exact flattened indices", {
  fixture <- make_tiny_crossed_design_fixture()
  truth <- fixture$truth

  expect_identical(unname(truth$state_id), rep(1:3, times = 4))
  expect_identical(unname(truth$stratum_id), rep(c(1L, 1L, 2L, 2L), each = 3))
  expect_identical(unname(truth$psu_flat_id), rep(1:4, each = 3))
  expect_identical(
    unname(truth$psu_within_stratum_id),
    rep(rep(1:2, each = 3), times = 2)
  )
  expect_identical(unname(truth$J_h), c(2L, 2L))
  expect_identical(unname(truth$psu_start), c(1L, 3L))
  expect_equal(sum(truth$J_h), truth$dimensions[["J"]])
  expect_identical(truth$composite_psu_map$flat_id, 1:4)
})

test_that("crossed-design weights, shares, and vhat have frozen truth", {
  fixture <- make_tiny_crossed_design_fixture()
  truth <- fixture$truth
  known <- fixture$known

  expect_identical(unname(truth$raw_weights), fixture$data$weight)
  expect_true(all(truth$raw_weights > 0))
  expect_gt(stats::sd(truth$raw_weights), 0)
  expect_equal(
    unname(truth$mean_one_weights),
    fixture$data$weight / mean(fixture$data$weight),
    tolerance = 1e-12
  )
  expect_equal(sum(truth$mean_one_weights), truth$dimensions[["N"]],
               tolerance = 1e-12)
  expect_equal(mean(truth$mean_one_weights), 1, tolerance = 1e-12)

  expect_named(known$population_shares, truth$state_levels)
  expect_equal(sum(known$population_shares), 1, tolerance = 1e-12)
  expect_true(all(known$population_shares > 0))
  expect_named(known$vhat, truth$state_levels)
  expect_true(all(is.finite(known$vhat)))
  expect_true(all(known$vhat >= 0))
})
