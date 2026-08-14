test_that("expanded G4-R truth-only DGP has frozen crossed dimensions", {
  truth <- generate_g4r_truth_only("low", 7310010L)
  expect_s3_class(truth, "bhf_g4r_truth_only")
  expect_identical(truth$truth$outcome_generated, FALSE)
  expect_false("outcome" %in% names(truth$data))
  expect_identical(
    unname(truth$truth$config$dimensions), c(960L, 20L, 8L, 24L)
  )
  expect_equal(length(unique(truth$data$state)), 20L)
  expect_equal(length(unique(truth$data$stratum)), 8L)
  expect_equal(length(unique(truth$data$psu)), 24L)
  expect_true(all(table(truth$truth$psu_structure$stratum_id) == 3L))
})

test_that("truth-only DGP is deterministic and leaves global RNG unchanged", {
  set.seed(875L)
  before <- .Random.seed
  first <- generate_g4r_truth_only("high", 7410010L)
  after <- .Random.seed
  second <- generate_g4r_truth_only("high", 7410010L)
  expect_identical(after, before)
  expect_identical(first, second)
  expect_identical(
    first$truth$checksums$data_without_outcome,
    second$truth$checksums$data_without_outcome
  )
})

test_that("expanded truth effects satisfy all centering invariants", {
  truth <- generate_g4r_truth_only("low", 7310020L)$truth
  expect_lt(abs(truth$centering_residuals$state_weighted), 1e-12)
  expect_lt(abs(truth$centering_residuals$stratum_mean), 1e-12)
  expect_true(all(abs(truth$centering_residuals$psu_within_stratum) < 1e-12))
  expect_equal(sum(truth$population_shares), 1, tolerance = 1e-14)
  expect_equal(truth$weights$likelihood_sum, 960, tolerance = 1e-10)
})

test_that("low and high profiles retain the published SD anchors", {
  low <- generate_g4r_truth_only("low", 7310030L)$truth
  high <- generate_g4r_truth_only("high", 7410030L)$truth
  expect_identical(low$sigmas, c(state=0.3,stratum=0.4,psu=0.5))
  expect_identical(high$sigmas, c(state=0.6,stratum=0.7,psu=0.9))
  expect_true(all(vapply(low$estimands[c("A","A_star","B")],
                         is.list, logical(1L))))
})

test_that("truth-only DGP rejects invalid vhat and never adds outcomes", {
  expect_error(
    generate_g4r_truth_only("low", 7310040L, vhat_state = rep(-1,20)),
    class = "bhf_g4r_preflight_error"
  )
  value <- generate_g4r_truth_only("low", 7310040L,
                                   vhat_state = rep(0.001,20))
  expect_false(any(grepl("outcome", names(value$data), fixed=TRUE)))
  expect_false("outcomes" %in% names(value$truth$checksums))
})

test_that("fast candidate gap oracle matches the full truth-only generator", {
  cases <- list(
    c(profile="low",seed=7310010L), c(profile="low",seed=7310020L),
    c(profile="high",seed=7410010L), c(profile="high",seed=7410020L)
  )
  for (case in cases) {
    profile <- case[["profile"]]
    seed <- as.integer(case[["seed"]])
    full <- generate_g4r_truth_only(profile,seed)$truth$gap
    fast <- bhf_g4r_candidate_gap_oracle(profile,seed)
    expect_lt(abs(fast$delta_truth-full$delta_truth),1e-12)
    expect_lt(abs(fast$r_truth-full$r_truth),1e-12)
    expect_identical(fast$regime,full$regime)
    expect_identical(fast$direction,full$direction)
  }
})
