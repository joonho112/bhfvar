test_that("prepared-data methods report scientific provenance", {
  fixture <- make_tiny_crossed_design_fixture()
  prepared <- suppressWarnings(prepare_bhf_data(
    fixture$data, "outcome", "state", "stratum", "psu", "weight",
    population_shares = fixture$known$population_shares,
    deattenuation = "none"
  ))
  printed <- capture.output(print(prepared))
  summarized <- capture.output(summary(prepared))
  expect_match(paste(printed, collapse = "\n"), "schema 0.5.0")
  expect_match(paste(printed, collapse = "\n"), "Likelihood weights: mean_one")
  expect_match(paste(printed, collapse = "\n"),
               "A*: none (unavailable)", fixed = TRUE)
  expect_match(paste(summarized, collapse = "\n"), "Scientific Inputs")
  expect_match(paste(summarized, collapse = "\n"),
               "half-t(3, 0, 2.5)", fixed = TRUE)
})

test_that("result methods expose estimand, interval, and limitation labels", {
  local_result_draw_backend()
  fit <- make_result_contract_fit(FALSE)
  decomposition <- variance_decomposition(fit, prob = 0.90)
  domains <- domain_estimates(fit, "B", prob = 0.90)
  overall <- overall_estimate(fit, "B", prob = 0.90)
  expect_warning(
    likelihood <- log_lik(fit, "pseudo", "psu"),
    class = "bhf_pseudo_loo_warning"
  )

  decomposition_text <- capture.output(print(decomposition))
  domain_text <- capture.output(print(domains))
  overall_text <- capture.output(print(overall))
  likelihood_text <- capture.output(print(likelihood))
  expect_match(paste(decomposition_text, collapse = "\n"),
               "A*: unavailable", fixed = TRUE)
  expect_match(paste(decomposition_text, collapse = "\n"),
               "90% credible interval")
  expect_match(paste(domain_text, collapse = "\n"), "Estimand B")
  expect_match(paste(overall_text, collapse = "\n"), "Estimand B")
  expect_match(paste(likelihood_text, collapse = "\n"),
               "not leave-one-unit refitting")
})

test_that("fit methods use versioned provenance and current result schema", {
  local_result_draw_backend()
  fit <- make_result_contract_fit(TRUE)
  fit$contract_id <- "bhfvar-fit-contract-0.5.0"
  fit$model_sha256 <- paste(rep("a", 64), collapse = "")
  fit$provenance <- list(sampling = list(
    chains = 2L, iter = 20L, warmup = 10L, seed = 42L,
    control = list(adapt_delta = 0.95, max_treedepth = 12L)
  ))
  fit$diagnostics <- list(
    n_divergent = 0L, max_depth_hits = 0L,
    rhat_max = 1.01, n_eff_min = 250
  )

  printed <- capture.output(print(fit))
  summarized <- capture.output(summary(fit))
  expect_match(paste(printed, collapse = "\n"), "schema 0.5.0")
  expect_match(paste(summarized, collapse = "\n"),
               "Article-aligned variance results")
  expect_match(paste(summarized, collapse = "\n"), "B_minus_A")
})
