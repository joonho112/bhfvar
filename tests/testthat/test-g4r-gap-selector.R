make_g4r_selector_ledger <- function() {
  rows <- list()
  index <- 0L
  for (profile in c("low", "high")) {
    for (candidate in seq_len(20L)) {
      index <- index + 1L
      regime <- if (candidate <= 4L) "grey" else if (candidate <= 10L) {
        "detectable"
      } else {
        "equivalence"
      }
      direction <- if (candidate <= 7L) "positive" else if (candidate <= 10L) {
        "negative"
      } else {
        "positive"
      }
      seed <- bhf_g4r_candidate_seeds(profile, candidate)
      rows[[index]] <- data.frame(
        profile = profile, candidate_index = candidate,
        effect_seed = unname(seed[["effect_seed"]]),
        delta_truth = if (direction == "negative") -0.2 else 0.2,
        r_truth = if (regime == "detectable") 0.2 else if (
          regime == "equivalence"
        ) 0.005 else 0.05,
        regime = regime, direction = direction,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

test_that("G4-R protocol freezes approved thresholds and dimensions", {
  protocol <- bhf_g4r_protocol()
  expect_identical(protocol$dgp$expected_n, 960L)
  expect_identical(protocol$gap$equivalence_truth_max, 0.01)
  expect_identical(protocol$gap$posterior_rope_max, 0.05)
  expect_identical(protocol$gap$detectable_truth_min, 0.10)
  expect_identical(protocol$rerun$eligibility,
                   "sampler_diagnostic_failure_only")
  expect_match(protocol$authority$protocol_sha256, "^[0-9a-f]{64}$")
})

test_that("gap metric separates equivalence, grey, and detectable truth", {
  make_estimands <- function(a, b) list(
    A = list(summary = c(between = a)),
    B = list(summary = c(between = b))
  )
  expect_identical(bhf_g4r_gap_metrics(make_estimands(1, 0.995))$regime,
                   "equivalence")
  expect_identical(bhf_g4r_gap_metrics(make_estimands(1, 0.95))$regime,
                   "grey")
  detectable <- bhf_g4r_gap_metrics(make_estimands(1, 0.8))
  expect_identical(detectable$regime, "detectable")
  expect_identical(detectable$direction, "positive")
  negative <- bhf_g4r_gap_metrics(make_estimands(1, 1.2))
  expect_identical(negative$direction, "negative")
})

test_that("candidate seeds are deterministic, prospective, and separated", {
  expect_identical(
    bhf_g4r_candidate_seeds("low", 1L),
    c(effect_seed = 7310010L, outcome_seed = 7310011L,
      mcmc_seed = 7310012L)
  )
  expect_identical(
    bhf_g4r_candidate_seeds("high", 5000L)[["effect_seed"]],
    7460000L
  )
  expect_error(bhf_g4r_candidate_seeds("low", 0L),
               class = "bhf_g4r_preflight_error")
})

test_that("selector uses first outcome-blind candidates in each quota", {
  ledger <- make_g4r_selector_ledger()
  selected <- bhf_g4r_select_candidates(ledger[rev(seq_len(nrow(ledger))), ])
  expect_equal(nrow(selected), 24L)
  counts <- table(selected$profile, selected$selection_group)
  expect_identical(unname(counts[, "detectable_positive"]), c(3L, 3L))
  expect_identical(unname(counts[, "detectable_negative"]), c(3L, 3L))
  expect_identical(unname(counts[, "equivalence"]), c(6L, 6L))
  low_positive <- subset(
    selected, profile == "low" & selection_group == "detectable_positive"
  )
  expect_identical(low_positive$candidate_index, 5:7)
  expect_false(any(c("outcome", "posterior", "draws") %in% names(selected)))
  expect_true(all(c("outcome_seed", "mcmc_seed") %in% names(selected)))
})

test_that("selector fails closed on outcome fields and unfilled quotas", {
  ledger <- make_g4r_selector_ledger()
  ledger$outcome <- 0L
  expect_error(bhf_g4r_select_candidates(ledger),
               class = "bhf_g4r_preflight_error")
  incomplete <- subset(make_g4r_selector_ledger(), candidate_index <= 12L)
  expect_error(bhf_g4r_select_candidates(incomplete),
               "did not fill", class = "bhf_g4r_preflight_error")
})
