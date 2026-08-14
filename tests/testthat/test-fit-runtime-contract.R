.bhf_runtime_prepared <- function() {
  fixture <- make_tiny_crossed_design_fixture()
  suppressWarnings(prepare_bhf_data(
    fixture$data,
    outcome = "outcome",
    domain = "state",
    strata = "stratum",
    psu = "psu",
    weights = "weight",
    population_shares = fixture$known$population_shares,
    deattenuation = "supplied",
    sampling_variances = fixture$known$vhat,
    sampling_variance_method = "external_taylor"
  ))
}

.bhf_runtime_model <- function(hash = paste(rep("a", 64L), collapse = "")) {
  model <- structure(list(fake = TRUE), class = "stanmodel")
  attr(model, "bhfvar_model_sha256") <- hash
  model
}

.bhf_runtime_summary_matrix <- function() {
  matrix(
    c(
      0.0, 0.1, 400, 1.00,
      0.5, 0.1, 350, 1.01,
      0.4, 0.1, 300, 1.00,
      0.3, 0.1, 250, 1.00,
      9.9, 0.1, 1, 9.99
    ),
    nrow = 5L,
    byrow = TRUE,
    dimnames = list(
      c("alpha", "sigma_state", "sigma_stratum", "sigma_psu",
        "var_total_A"),
      c("mean", "sd", "n_eff", "Rhat")
    )
  )
}

test_that("bhf_fit returns a versioned, source-fingerprinted provenance object", {
  prepared <- .bhf_runtime_prepared()
  model <- .bhf_runtime_model()
  state <- new.env(parent = emptyenv())
  state$auto_write <- FALSE
  state$sampling_args <- NULL
  state$summary_pars <- NULL

  old_mc_present <- "mc.cores" %in% names(options())
  old_mc_value <- getOption("mc.cores")
  on.exit({
    if (old_mc_present) options(mc.cores = old_mc_value) else options(mc.cores = NULL)
  }, add = TRUE)
  options(mc.cores = 7L)

  local_mocked_bindings(
    .bhf_get_rstan_auto_write = function() state$auto_write,
    .bhf_set_rstan_auto_write = function(value) {
      state$auto_write <- value
      invisible(value)
    },
    .bhf_rstan_sampling = function(...) {
      state$sampling_args <- list(...)
      state$auto_write_during_sampling <- state$auto_write
      state$mc_cores_during_sampling <- getOption("mc.cores")
      structure(list(fake = TRUE), class = "stanfit")
    },
    .bhf_rstan_sampler_params = function(...) {
      list(
        cbind(divergent__ = c(0, 0, 0), treedepth__ = c(5, 6, 8)),
        cbind(divergent__ = c(0, 0, 0), treedepth__ = c(6, 7, 8))
      )
    },
    .bhf_rstan_summary = function(object, pars, ...) {
      state$summary_pars <- pars
      list(summary = .bhf_runtime_summary_matrix())
    },
    .package = "bhfvar"
  )

  fit <- suppressMessages(bhf_fit(
    prepared,
    model = model,
    chains = 2,
    iter = 20,
    warmup = 10,
    seed = 90210,
    cores = 2,
    adapt_delta = 0.93,
    max_treedepth = 9,
    refresh = 0,
    algorithm = "NUTS"
  ))

  expect_s3_class(fit, "bhf_fit")
  expect_identical(fit$schema_version, "0.5.0")
  expect_identical(fit$contract_id, "bhfvar-fit-contract-0.5.0")
  expect_true(is.character(fit$package_version) &&
                length(fit$package_version) == 1L)
  expect_match(fit$model_sha256, "^[0-9a-f]{64}$")
  expect_identical(fit$model_sha256,
                   attr(model, "bhfvar_model_sha256", exact = TRUE))
  expect_identical(fit$provenance$contract$schema_version, "0.5.0")
  expect_identical(fit$provenance$contract$data_contract_id,
                   prepared$contract_id)
  expect_identical(fit$provenance$contract$model_sha256,
                   fit$model_sha256)
  expect_identical(fit$provenance$data$dimensions,
                   prepared$stan_data[c("N", "S", "H", "J")])
  expect_identical(fit$provenance$weights$raw, prepared$weight_info$raw)
  expect_identical(fit$provenance$weights$w_lik,
                   prepared$weight_info$likelihood)
  expect_identical(fit$provenance$weights$method,
                   prepared$provenance$weights$method)
  expect_identical(fit$provenance$population_shares$source,
                   prepared$provenance$population_shares$source)
  expect_identical(fit$provenance$population_shares$values,
                   prepared$population_share_info$values)
  expect_identical(fit$provenance$sampling_variances$named_values,
                   prepared$sampling_variance_info$named_values)
  expect_identical(fit$provenance$sampling_variances$stan_values,
                   prepared$sampling_variance_info$stan_values)
  expect_identical(fit$provenance$prior, prepared$provenance$prior)

  sampling <- fit$provenance$sampling
  expect_identical(sampling$seed, 90210L)
  expect_identical(sampling$control,
                   list(adapt_delta = 0.93, max_treedepth = 9L))
  expect_identical(sampling$additional_argument_names, "algorithm")
  expect_identical(sampling$additional_arguments, list(algorithm = "NUTS"))
  expect_identical(state$sampling_args$control, sampling$control)
  expect_identical(state$sampling_args$algorithm, "NUTS")
  expect_true(state$auto_write_during_sampling)
  expect_identical(state$mc_cores_during_sampling, 2L)
  expect_false(state$auto_write)
  expect_identical(getOption("mc.cores"), 7L)

  session <- fit$provenance$session
  expect_true(all(c(
    "started_at_utc", "completed_at_utc", "r_version", "platform", "os",
    "release", "bhfvar_version", "rstan_version"
  ) %in% names(session)))
  expect_true(all(nzchar(unlist(session[c(
    "started_at_utc", "completed_at_utc", "r_version", "platform"
  )]))))
  expect_identical(
    state$summary_pars,
    c("alpha", "sigma_state", "sigma_stratum", "sigma_psu")
  )
  expect_identical(fit$diagnostics$parameters,
                   c("alpha", "sigma_state", "sigma_stratum", "sigma_psu"))
  expect_false(any(grepl("ess_bulk|ess_tail", names(fit$diagnostics))))
  expect_identical(
    fit$diagnostics$n_eff_label,
    "legacy rstan n_eff (not bulk or tail ESS)"
  )
})

test_that("sampling failure restores auto_write and an absent mc.cores option", {
  prepared <- .bhf_runtime_prepared()
  model <- .bhf_runtime_model()
  state <- new.env(parent = emptyenv())
  state$auto_write <- FALSE

  old_mc_present <- "mc.cores" %in% names(options())
  old_mc_value <- getOption("mc.cores")
  on.exit({
    if (old_mc_present) options(mc.cores = old_mc_value) else options(mc.cores = NULL)
  }, add = TRUE)
  options(mc.cores = NULL)

  local_mocked_bindings(
    .bhf_get_rstan_auto_write = function() state$auto_write,
    .bhf_set_rstan_auto_write = function(value) {
      state$auto_write <- value
      invisible(value)
    },
    .bhf_rstan_sampling = function(...) stop("forced sampling failure"),
    .package = "bhfvar"
  )

  expect_error(
    suppressMessages(bhf_fit(
      prepared,
      model = model,
      chains = 1,
      iter = 10,
      warmup = 5,
      cores = 1,
      refresh = 0
    )),
    "forced sampling failure"
  )
  expect_false(state$auto_write)
  expect_false("mc.cores" %in% names(options()))
})

test_that("compile success and failure both restore rstan auto_write", {
  stan_candidates <- c(
    testthat::test_path("..", "..", "inst", "stan", "bhf_hybrid.stan"),
    system.file("stan", "bhf_hybrid.stan", package = "bhfvar")
  )
  stan_file <- normalizePath(
    stan_candidates[file.exists(stan_candidates)][[1L]], mustWork = TRUE
  )

  success_state <- new.env(parent = emptyenv())
  success_state$auto_write <- FALSE
  local_mocked_bindings(
    .bhf_get_rstan_auto_write = function() success_state$auto_write,
    .bhf_set_rstan_auto_write = function(value) {
      success_state$auto_write <- value
      invisible(value)
    },
    .bhf_locate_stan_file = function() stan_file,
    .bhf_rstan_stan_model = function(...) {
      success_state$during_compile <- success_state$auto_write
      structure(list(fake = TRUE), class = "stanmodel")
    },
    .package = "bhfvar"
  )
  model <- compile_bhf_model(verbose = FALSE, auto_write = TRUE)
  expect_true(success_state$during_compile)
  expect_false(success_state$auto_write)
  expect_match(attr(model, "bhfvar_model_sha256", exact = TRUE),
               "^[0-9a-f]{64}$")
  expect_identical(attr(model, "bhfvar_model_path", exact = TRUE), stan_file)

  failure_state <- new.env(parent = emptyenv())
  failure_state$auto_write <- TRUE
  local_mocked_bindings(
    .bhf_get_rstan_auto_write = function() failure_state$auto_write,
    .bhf_set_rstan_auto_write = function(value) {
      failure_state$auto_write <- value
      invisible(value)
    },
    .bhf_locate_stan_file = function() stan_file,
    .bhf_rstan_stan_model = function(...) stop("forced compile failure"),
    .package = "bhfvar"
  )
  expect_error(
    compile_bhf_model(verbose = FALSE, auto_write = FALSE),
    "Stan model compilation failed.*forced compile failure"
  )
  expect_true(failure_state$auto_write)
})

test_that("fit and compile scalar arguments fail before runtime mutation", {
  prepared <- .bhf_runtime_prepared()
  model <- .bhf_runtime_model()
  invalid <- list(
    list(chains = 0),
    list(chains = 1.5),
    list(iter = c(10, 20)),
    list(iter = 0),
    list(iter = 10, warmup = 10),
    list(seed = 0),
    list(seed = NA_real_),
    list(cores = 0),
    list(adapt_delta = 0),
    list(adapt_delta = 1),
    list(adapt_delta = NA_real_),
    list(max_treedepth = 2.5),
    list(refresh = -1)
  )
  for (arguments in invalid) {
    expect_error(
      suppressMessages(do.call(
        bhf_fit,
        c(list(data = prepared, model = model), arguments)
      )),
      info = paste(names(arguments), collapse = ",")
    )
  }

  expect_error(
    suppressMessages(bhf_fit(
      prepared, model = model, control = list(adapt_delta = 0.99)
    )),
    "reserved sampling argument.*control"
  )
  duplicated_dots <- structure(list(1, 2), names = c("foo", "foo"))
  expect_error(
    suppressMessages(do.call(
      bhf_fit,
      c(list(data = prepared, model = model), duplicated_dots)
    )),
    "unique names"
  )
  expect_error(compile_bhf_model(verbose = NA), "verbose")
  expect_error(compile_bhf_model(auto_write = c(TRUE, FALSE)), "auto_write")
})

test_that("diagnostics use configured treedepth and only structural parameters", {
  state <- new.env(parent = emptyenv())
  local_mocked_bindings(
    .bhf_rstan_sampler_params = function(...) {
      list(
        cbind(divergent__ = c(0, 1, 0), treedepth__ = c(8, 9, 10)),
        cbind(divergent__ = c(1, 0), treedepth__ = c(7, 9))
      )
    },
    .bhf_rstan_summary = function(object, pars, ...) {
      state$pars <- pars
      summary <- .bhf_runtime_summary_matrix()
      summary["sigma_stratum", "Rhat"] <- Inf
      summary["sigma_psu", "n_eff"] <- NA_real_
      list(summary = summary)
    },
    .package = "bhfvar"
  )
  pars <- c("alpha", "sigma_state", "sigma_stratum", "sigma_psu")
  diagnostics <- compute_diagnostics(
    structure(list(), class = "stanfit"),
    max_treedepth = 9,
    pars = pars
  )

  expect_identical(state$pars, pars)
  expect_identical(diagnostics$parameters, pars)
  expect_identical(diagnostics$n_divergent, 2L)
  expect_identical(diagnostics$max_treedepth_configured, 9L)
  expect_identical(diagnostics$max_treedepth_observed, 10L)
  expect_identical(diagnostics$max_depth_hits, 3L)
  expect_equal(diagnostics$rhat_max, 1.01)
  expect_identical(diagnostics$rhat_n, 3L)
  expect_equal(diagnostics$n_eff_min, 300)
  expect_identical(diagnostics$n_eff_n, 3L)
  expect_false(any(c("var_total_A", "p_bar_A") %in%
                     diagnostics$parameters))
})

test_that("empty diagnostics are missing-safe and never infinite", {
  local_mocked_bindings(
    .bhf_rstan_sampler_params = function(...) list(matrix(numeric(), 0, 0)),
    .bhf_rstan_summary = function(...) {
      list(summary = matrix(numeric(), 0, 0))
    },
    .package = "bhfvar"
  )
  diagnostics <- compute_diagnostics(
    structure(list(), class = "stanfit"),
    max_treedepth = 11
  )

  expect_identical(diagnostics$parameters, NULL)
  expect_identical(diagnostics$n_divergent, 0L)
  expect_identical(diagnostics$max_depth_hits, 0L)
  expect_true(is.na(diagnostics$max_treedepth_observed))
  expect_true(is.na(diagnostics$rhat_max))
  expect_true(is.na(diagnostics$rhat_all_ok))
  expect_true(is.na(diagnostics$n_eff_min))
  expect_true(is.na(diagnostics$n_eff_median))
  numeric_diagnostics <- unlist(diagnostics[vapply(
    diagnostics,
    function(x) is.numeric(x) && length(x) == 1L,
    logical(1)
  )])
  expect_false(any(is.infinite(numeric_diagnostics)))
  expect_false(any(is.nan(numeric_diagnostics)))
})
