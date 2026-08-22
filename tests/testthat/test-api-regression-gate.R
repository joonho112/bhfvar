test_that("0.5.0 public signatures and defaults are frozen", {
  expect_identical(
    names(formals(prepare_bhf_data)),
    c(
      "data", "outcome", "domain", "strata", "psu", "weights",
      "population_shares", "weight_scaling", "deattenuation",
      "sampling_variances", "sampling_variance_method",
      "use_deattenuation", "prior_alpha_mean", "prior_alpha_sd",
      "sigma_state_prior"
    )
  )
  expect_identical(
    names(formals(bhf_fit)),
    c(
      "data", "model", "chains", "iter", "warmup", "seed", "cores",
      "adapt_delta", "max_treedepth", "refresh", "..."
    )
  )
  expect_identical(eval(formals(bhf_fit)$adapt_delta), 0.95)
  expect_identical(eval(formals(variance_decomposition)$prob), 0.95)
  expect_identical(eval(formals(domain_estimates)$prob), 0.95)
  expect_identical(eval(formals(overall_estimate)$prob), 0.95)
  expect_identical(
    names(formals(log_lik)), c("fit", "kind", "aggregate")
  )
})

test_that("public result APIs expose only article-aligned schemas", {
  local_result_draw_backend()
  fit <- make_result_contract_fit(TRUE)
  decomposition <- variance_decomposition(fit)
  domain_a <- domain_estimates(fit, "A")
  domain_b <- domain_estimates(fit, "B")
  overall_a <- overall_estimate(fit, "A")
  overall_b <- overall_estimate(fit, "B")

  expect_named(
    decomposition,
    c("schema_version", "interval", "latent", "A", "A_star", "B",
      "gaps", "summary_table")
  )
  expect_false(any(c("logit", "prob", "deatten") %in%
                     names(decomposition)))
  expect_true(decomposition$A_star$available)
  expect_s3_class(domain_a, "bhf_domain_estimates")
  expect_s3_class(domain_b, "bhf_domain_estimates")
  expect_identical(unique(domain_a$estimand), "A")
  expect_identical(unique(domain_b$estimand), "B")
  expect_s3_class(overall_a, "bhf_overall_estimate")
  expect_s3_class(overall_b, "bhf_overall_estimate")
  expect_identical(overall_a$estimand, "A")
  expect_identical(overall_b$estimand, "B")
})

test_that("legacy fixtures fail closed through public fit and result paths", {
  legacy_data <- dget(testthat::test_path(
    "..", "fixtures", "bhf_data-0.3.0.dput"
  ))
  legacy_fit <- dget(testthat::test_path(
    "..", "fixtures", "bhf_fit-0.3.0.dput"
  ))
  fake_model <- structure(list(), class = "stanmodel")

  expect_error(
    bhf_fit(legacy_data, model = fake_model, chains = 1, iter = 2,
            warmup = 1, cores = 1, refresh = 0),
    class = "bhf_legacy_object_error"
  )
  for (extractor in list(
    variance_decomposition,
    function(x) domain_estimates(x, "A"),
    function(x) overall_estimate(x, "A"),
    function(x) log_lik(x, "raw")
  )) {
    error <- tryCatch(extractor(legacy_fit), error = identity)
    expect_s3_class(error, "bhf_legacy_object_error")
    expect_identical(error$detected_schema, "0.3.0")
    expect_identical(error$required_action, "refit")
  }
  data_print_error <- tryCatch(print(legacy_data), error = identity)
  expect_s3_class(data_print_error, "bhf_legacy_object_error")
  expect_identical(data_print_error$detected_schema, "0.3.0")
  data_summary_error <- tryCatch(summary(legacy_data), error = identity)
  expect_s3_class(data_summary_error, "bhf_legacy_object_error")
  fit_print_error <- tryCatch(print(legacy_fit), error = identity)
  expect_s3_class(fit_print_error, "bhf_legacy_object_error")
})

test_that("namespace surface registers all current public methods", {
  namespace_path <- testthat::test_path("..", "..", "NAMESPACE")
  exports <- c(
    "bhf_fit", "calc_eff_n", "compile_bhf_model", "domain_estimates",
    "get_stan_file_path", "log_lik", "overall_estimate",
    "prepare_bhf_data", "variance_decomposition",
    "bhf_plot_variance", "bhf_plot_icc", "bhf_plot_domains",
    "bhf_plot_shrinkage", "bhf_plot_astar_sensitivity",
    "bhf_astar_sensitivity"
  )
  if (!file.exists(namespace_path)) {
    expect_setequal(getNamespaceExports("bhfvar"), exports)
    registered <- getNamespaceInfo(asNamespace("bhfvar"), "S3methods")
    expect_true(all(c(
      "bhf_data", "bhf_fit", "bhf_variance_decomposition",
      "bhf_domain_estimates", "bhf_overall_estimate", "bhf_log_lik"
    ) %in% registered[, 2L]))
    return(invisible(NULL))
  }
  namespace <- readLines(namespace_path, warn = FALSE)
  for (name in exports) {
    expect_true(paste0("export(", name, ")") %in% namespace, info = name)
  }
  methods <- c(
    "bhf_data", "bhf_fit", "bhf_variance_decomposition",
    "bhf_domain_estimates", "bhf_overall_estimate", "bhf_log_lik"
  )
  for (class in methods) {
    expect_true(
      paste0("S3method(print,", class, ")") %in% namespace,
      info = class
    )
  }
})

test_that("legacy scientific names occur only in explicit migration inventory", {
  r_dir <- testthat::test_path("..", "..", "R")
  if (!dir.exists(r_dir)) {
    succeed("R sources are inspected before installed-test isolation")
    return(invisible(NULL))
  }
  files <- setdiff(list.files(r_dir, full.names = TRUE, pattern = "\\.R$"),
                   file.path(r_dir, "legacy.R"))
  source <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")
  forbidden <- c(
    "var_between_state", "p_state_marginal", "p_state_conditional",
    "icc_prob", "icc_deatten", "reliability_state", "c_factor"
  )
  for (name in forbidden) {
    expect_false(grepl(paste0("\\b", name, "\\b"), source, perl = TRUE),
                 info = name)
  }
})
