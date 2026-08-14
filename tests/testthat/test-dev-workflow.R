.bhf_dev_workflow_files <- function() {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  file.path(
    root, "dev",
    c(
      "00_verify_environment.R", "01_compile_model.R",
      "02_prepare_data.R", "03_fit_model.R", "04_extract_results.R",
      "05_install_and_test_package.R", "run_workflow.R"
    )
  )
}

.bhf_dev_with_environment <- function(values, code) {
  names <- names(values)
  old <- Sys.getenv(names, unset = NA_character_, names = TRUE)
  on.exit({
    for (name in names) {
      if (is.na(old[[name]])) {
        Sys.unsetenv(name)
      } else {
        do.call(Sys.setenv, stats::setNames(list(old[[name]]), name))
      }
    }
  }, add = TRUE)
  do.call(Sys.setenv, as.list(values))
  force(code)
}

.bhf_dev_sources_available <- function() {
  all(file.exists(.bhf_dev_workflow_files()))
}

.bhf_dev_archive_pass <- function() {
  succeed("dev/ is intentionally excluded from the source archive")
  invisible(NULL)
}

test_that("all dev workflow scripts parse", {
  if (!.bhf_dev_sources_available()) return(.bhf_dev_archive_pass())
  files <- .bhf_dev_workflow_files()
  expect_length(files, 7L)
  expect_true(all(file.exists(files)))
  for (path in files) expect_no_error(parse(file = path))
})

test_that("dev workflow contains current public contracts and no stale names", {
  if (!.bhf_dev_sources_available()) return(.bhf_dev_archive_pass())
  files <- .bhf_dev_workflow_files()
  source <- paste(
    unlist(lapply(files, readLines, warn = FALSE), use.names = FALSE),
    collapse = "\n"
  )

  required <- c(
    "schema_version", "weight_scaling = \"mean_one\"",
    "sigma_stratum", "sigma_state_prior = \"half_t3_2.5\"",
    "variance_decomposition", "estimand = \"A\"", "estimand = \"B\"",
    "kind = \"pseudo\"", "kind = \"raw\"",
    "ordinary LOO/WAIC", "provenance"
  )
  for (text in required) expect_match(source, text, fixed = TRUE, info = text)

  forbidden <- c(
    "(?i)2-level", "\\bc_factor\\b", "\\bp_state_marginal\\b",
    "\\breliability(_state|_avg)?\\b", "\\blegacy_d2\\b",
    "\\buse_deattenuation\\b", "\\bicc_prob\\b",
    "\\bicc_deatten\\b", "\\bvar_between_state\\b"
  )
  for (pattern in forbidden) {
    expect_false(grepl(pattern, source, perl = TRUE), info = pattern)
  }
})

test_that("heavy operations and artifact writes require explicit opt-in", {
  if (!.bhf_dev_sources_available()) return(.bhf_dev_archive_pass())
  source <- paste(
    unlist(lapply(.bhf_dev_workflow_files(), readLines, warn = FALSE),
           use.names = FALSE),
    collapse = "\n"
  )
  flags <- c(
    "BHFVAR_DEV_COMPILE", "BHFVAR_DEV_FIT", "BHFVAR_DEV_EXTRACT",
    "BHFVAR_DEV_WRITE", "BHFVAR_DEV_INSTALL_TEST"
  )
  for (flag in flags) expect_match(source, flag, fixed = TRUE, info = flag)
  expect_match(
    source,
    'if (.bhf_dev_flag("BHFVAR_DEV_COMPILE"))',
    fixed = TRUE
  )
  expect_match(
    source,
    'if (.bhf_dev_flag("BHFVAR_DEV_FIT"))',
    fixed = TRUE
  )
  expect_match(
    source,
    'if (.bhf_dev_flag("BHFVAR_DEV_INSTALL_TEST"))',
    fixed = TRUE
  )
})

test_that("default workflow is light, current, and does not refresh artifacts", {
  if (!.bhf_dev_sources_available()) return(.bhf_dev_archive_pass())
  files <- .bhf_dev_workflow_files()
  root <- dirname(dirname(files[[1L]]))
  artifact_paths <- list.files(
    file.path(root, "dev", "data"), full.names = TRUE, recursive = TRUE
  )
  before_hash <- if (length(artifact_paths)) tools::md5sum(artifact_paths) else
    character()
  before_time <- if (length(artifact_paths)) file.info(artifact_paths)$mtime else
    as.POSIXct(character())
  working_directory <- getwd()

  values <- c(
    BHFVAR_DEV_ROOT = root,
    BHFVAR_DEV_COMPILE = "false",
    BHFVAR_DEV_FIT = "false",
    BHFVAR_DEV_EXTRACT = "false",
    BHFVAR_DEV_WRITE = "false",
    BHFVAR_DEV_INSTALL_TEST = "false",
    BHFVAR_DEV_DEATTENUATION = "none"
  )
  workflow_environment <- new.env(parent = globalenv())
  output <- .bhf_dev_with_environment(values, {
    capture.output(
      expect_no_error(sys.source(files[[7L]], envir = workflow_environment))
    )
  })

  expect_true(exists("prepared_data", envir = workflow_environment,
                     inherits = FALSE))
  prepared <- get("prepared_data", envir = workflow_environment,
                  inherits = FALSE)
  expect_s3_class(prepared, "bhf_data")
  expect_identical(prepared$schema_version, "0.5.0")
  expect_identical(prepared$provenance$weights$method, "mean_one")
  expect_equal(sum(prepared$stan_data$w_lik), prepared$stan_data$N,
               tolerance = 1e-12)
  expect_true("psu_flat_id" %in% names(prepared$stan_data))
  expect_identical(prepared$prior_info$sigma_state$variant, "half_t3_2.5")
  expect_identical(prepared$provenance$sampling_variances$mode, "none")
  expect_false(exists("bhf_model", envir = workflow_environment,
                      inherits = FALSE))
  expect_false(exists("bhf_fit_result", envir = workflow_environment,
                      inherits = FALSE))
  expect_false(exists("bhf_results", envir = workflow_environment,
                      inherits = FALSE))
  expect_true(any(grepl("No generated result or model artifact was written",
                        output, fixed = TRUE)))
  expect_identical(getwd(), working_directory)

  after_paths <- list.files(
    file.path(root, "dev", "data"), full.names = TRUE, recursive = TRUE
  )
  expect_identical(after_paths, artifact_paths)
  if (length(artifact_paths)) {
    expect_identical(unname(tools::md5sum(after_paths)), unname(before_hash))
    expect_identical(file.info(after_paths)$mtime, before_time)
  }
})

test_that("environment bootstrap loads every workflow public API", {
  if (!.bhf_dev_sources_available()) return(.bhf_dev_archive_pass())
  files <- .bhf_dev_workflow_files()
  root <- dirname(dirname(files[[1L]]))
  environment <- new.env(parent = globalenv())
  output <- .bhf_dev_with_environment(c(BHFVAR_DEV_ROOT = root), {
    capture.output(expect_no_error(sys.source(files[[1L]],
                                               envir = environment)))
  })
  expected <- c(
    "prepare_bhf_data", "compile_bhf_model", "bhf_fit",
    "variance_decomposition", "domain_estimates", "overall_estimate",
    "log_lik", "detect_bhf_object_schema"
  )
  expect_true(all(vapply(
    expected, exists, logical(1), envir = environment, mode = "function",
    inherits = FALSE
  )))
  expect_true(any(grepl("prepared-data schema: 0.5.0", output, fixed = TRUE)))
})
