test_that("G4-R weight arm produces a hard nondegenerate contrast", {
  truth <- generate_g4r_truth_with_vhat("low",7310010L)
  arm <- bhf_g4r_weight_arm(truth)
  expect_identical(arm$metrics$contrast_pass,TRUE)
  expect_gte(arm$metrics$max_abs_delta,0.10)
  expect_gte(arm$metrics$mean_abs_delta,0.02)
  expect_equal(arm$metrics$mean_one_sum,960,tolerance=1e-10)
  expect_equal(arm$metrics$legacy_d2_sum,960,tolerance=1e-10)
  expect_true(all(arm$raw_weight>0 & is.finite(arm$raw_weight)))
  expect_true(any(grepl("deprecated",arm$warning_messages)))
})

test_that("weight arm is deterministic and row-permutation equivariant", {
  truth <- generate_g4r_truth_with_vhat("high",7410020L)
  first <- bhf_g4r_weight_arm(truth)
  second <- bhf_g4r_weight_arm(truth)
  expect_identical(first,second)
  order <- rev(seq_len(nrow(truth$data)))
  permuted <- truth
  permuted$data <- permuted$data[order,,drop=FALSE]
  observed <- bhf_g4r_weight_arm(permuted)
  expect_equal(observed$raw_weight[order(order)],first$raw_weight,tolerance=0)
  expect_equal(observed$mean_one[order(order)],first$mean_one,tolerance=1e-14)
  expect_equal(observed$legacy_d2[order(order)],first$legacy_d2,tolerance=1e-14)
})

test_that("both scaling paths are invariant to global raw-weight scale", {
  truth <- generate_g4r_truth_with_vhat("low",7310030L)
  arm <- bhf_g4r_weight_arm(truth)
  state_id <- truth$data$state_id
  scaled_mean <- scale_likelihood_weights(100*arm$raw_weight,"mean_one",state_id)
  warnings <- character()
  scaled_legacy <- withCallingHandlers(
    scale_likelihood_weights(100*arm$raw_weight,"legacy_d2",state_id),
    warning=function(c) { warnings <<- c(warnings,conditionMessage(c));
                          invokeRestart("muffleWarning") }
  )
  expect_equal(scaled_mean,arm$mean_one,tolerance=1e-12)
  expect_equal(scaled_legacy,arm$legacy_d2,tolerance=1e-12)
  expect_true(any(grepl("deprecated",warnings)))
})

test_that("weight arm changes within-state profiles and rejects bad input", {
  truth <- generate_g4r_truth_with_vhat("low",7310040L)
  arm <- bhf_g4r_weight_arm(truth)
  cv <- vapply(split(arm$raw_weight,truth$data$state_id),function(x) {
    stats::sd(x)/mean(x)
  },numeric(1L))
  expect_gt(length(unique(round(cv,8))),1L)
  expect_error(bhf_g4r_weight_arm(list()),
               class="bhf_g4r_preflight_error")
})
