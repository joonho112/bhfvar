g4r_direct_astar_oracle <- function(truth) {
  shares <- truth$population_shares
  p <- stats::plogis(truth$config$alpha + truth$effects$state)
  mean_p <- sum(shares * p)
  between <- sum(shares * (p - mean_p)^2)
  within <- sum(shares * p * (1-p))
  correction <- sum(shares * truth$vhat_state)
  between_star <- max(0, between - correction)
  total_star <- between_star + within
  c(
    mean=mean_p, between=between_star, within=within, total=total_star,
    proportion=between_star/total_star, correction=correction
  )
}

test_that("nonzero vhat has the frozen heterogeneous weighted correction", {
  result <- generate_g4r_truth_with_vhat("low", 7310010L)
  truth <- result$truth
  design <- truth$vhat_design
  expect_true(all(truth$vhat_state > 0))
  expect_gt(length(unique(round(truth$vhat_state, 15))), 1L)
  expect_identical(names(truth$vhat_state), names(truth$population_shares))
  expect_equal(design$observed_correction, design$target_correction,
               tolerance=1e-14)
  expect_equal(
    design$target_correction,
    0.25 * truth$estimands$A$summary[["between"]], tolerance=1e-14
  )
})

test_that("A-star truth is independently 75 percent of A between", {
  for (profile in c("low","high")) {
    base <- if (profile == "low") 7310000L else 7410000L
    truth <- generate_g4r_truth_with_vhat(profile, base+20L)$truth
    oracle <- g4r_direct_astar_oracle(truth)
    expect_equal(
      unname(truth$estimands$A_star$summary[c(
        "mean","between","within","total","proportion"
      )]),
      unname(oracle[c("mean","between","within","total","proportion")]),
      tolerance=1e-12
    )
    expect_equal(
      truth$estimands$A_star$summary[["between"]],
      0.75 * truth$estimands$A$summary[["between"]], tolerance=1e-12
    )
    expect_false(truth$estimands$A_star$at_boundary)
    expect_false(truth$estimands$A_star$truncated)
  }
})

test_that("adding vhat changes A-star truth but no DGP row or A/B truth", {
  zero <- generate_g4r_truth_only("low",7310030L)
  nonzero <- generate_g4r_truth_with_vhat("low",7310030L)
  expect_identical(zero$data,nonzero$data)
  expect_identical(zero$truth$estimands$A,nonzero$truth$estimands$A)
  expect_identical(zero$truth$estimands$B,nonzero$truth$estimands$B)
  expect_false(identical(zero$truth$estimands$A_star,
                         nonzero$truth$estimands$A_star))
  expect_false("outcome" %in% names(nonzero$data))
})

test_that("vhat constructor fails closed for incomplete or degenerate truth", {
  expect_error(bhf_g4r_vhat_from_truth(list()),
               class="bhf_g4r_preflight_error")
  bad <- generate_g4r_truth_only("low",7310040L)$truth
  bad$estimands$A$summary[["between"]] <- 0
  expect_error(bhf_g4r_vhat_from_truth(bad),
               class="bhf_g4r_preflight_error")
})
