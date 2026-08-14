test_that("nested PSU index matches the crossed-design fixture truth", {
  fixture <- make_tiny_crossed_design_fixture()
  truth <- fixture$truth

  index <- build_nested_psu_index(
    unname(truth$stratum_id),
    fixture$data$psu
  )

  expect_named(index, c("J", "J_h", "start", "local", "flat", "mapping"))
  expect_identical(index$J, truth$dimensions[["J"]])
  expect_identical(index$J_h, unname(truth$J_h))
  expect_identical(index$start, unname(truth$psu_start))
  expect_identical(index$local, unname(truth$psu_within_stratum_id))
  expect_identical(index$flat, unname(truth$psu_flat_id))
  expect_identical(
    index$mapping,
    data.frame(
      stratum_id = c(1L, 1L, 2L, 2L),
      psu_label = c("P1", "P2", "P1", "P2"),
      local = c(1L, 2L, 1L, 2L),
      flat = 1:4,
      stringsAsFactors = FALSE
    )
  )
})

test_that("nested PSU index is invariant to row permutation", {
  fixture <- make_tiny_crossed_design_fixture()
  stratum_id <- unname(fixture$truth$stratum_id)
  psu <- fixture$data$psu
  reference <- build_nested_psu_index(stratum_id, psu)
  permutation <- c(12L, 1L, 7L, 4L, 10L, 2L, 8L, 5L, 11L, 3L, 9L, 6L)

  permuted <- build_nested_psu_index(
    stratum_id[permutation],
    psu[permutation]
  )

  expect_identical(permuted$J, reference$J)
  expect_identical(permuted$J_h, reference$J_h)
  expect_identical(permuted$start, reference$start)
  expect_identical(permuted$mapping, reference$mapping)
  expect_identical(permuted$local[order(permutation)], reference$local)
  expect_identical(permuted$flat[order(permutation)], reference$flat)
})

test_that("nested PSU index handles separator-collision labels without pasting", {
  collision_pairs <- data.frame(
    stratum = c("A::B", "A"),
    psu = c("C", "B::C"),
    stringsAsFactors = FALSE
  )
  stratum_levels <- sort(unique(collision_pairs$stratum))
  stratum_id <- match(collision_pairs$stratum, stratum_levels)

  expect_identical(
    paste(collision_pairs$stratum, collision_pairs$psu, sep = "::"),
    rep("A::B::C", 2)
  )

  index <- build_nested_psu_index(stratum_id, collision_pairs$psu)
  expect_identical(index$J, 2L)
  expect_identical(index$J_h, c(1L, 1L))
  expect_identical(sort(index$flat), c(1L, 2L))
  expect_identical(index$mapping$psu_label, c("B::C", "C"))
})

test_that("repeated observations share local and flat PSU indices", {
  stratum_id <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)
  psu <- c("P1", "P1", "P2", "P2", "P1", "P1", "P2", "P2")

  index <- build_nested_psu_index(stratum_id, psu)

  expect_identical(index$local, c(1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L))
  expect_identical(index$flat, c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L))
  expect_identical(index$J, 4L)
})

test_that("prepare_bhf_data preserves nested local and flat PSU indices", {
  fixture <- make_tiny_crossed_design_fixture()

  expect_warning(
    prepared <- prepare_bhf_data(
      fixture$data,
      outcome = "outcome",
      domain = "state",
      strata = "stratum",
      psu = "psu",
      weights = "weight",
      population_shares = fixture$known$population_shares
    ),
    "Sample size \\(12\\) is very small"
  )

  expect_identical(prepared$stan_data$J, fixture$truth$dimensions[["J"]])
  expect_identical(as.vector(prepared$stan_data$J_h),
                   unname(fixture$truth$J_h))
  expect_identical(as.vector(prepared$stan_data$psu_start),
                   unname(fixture$truth$psu_start))
  expect_identical(as.vector(prepared$stan_data$psu_in_stratum_id),
                   unname(fixture$truth$psu_within_stratum_id))
  expect_identical(as.vector(prepared$stan_data$psu_flat_id),
                   unname(fixture$truth$psu_flat_id))
  expect_identical(prepared$mapping$psu$label, c("P1", "P2", "P1", "P2"))
  expect_identical(prepared$mapping$psu$id, 1:4)
  expect_identical(prepared$mapping$psu$stratum_id, c(1L, 1L, 2L, 2L))
  expect_identical(prepared$mapping$psu$within_stratum_id,
                   c(1L, 2L, 1L, 2L))
})
