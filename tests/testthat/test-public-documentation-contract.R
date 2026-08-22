read_public_text <- function() {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  paths <- c(
    file.path(root, "DESCRIPTION"),
    file.path(root, "README.md"),
    file.path(root, "NEWS.md"),
    file.path(root, "R", "bhfvar-package.R"),
    list.files(file.path(root, "vignettes"), pattern = "[.]Rmd$",
               full.names = TRUE),
    file.path(root, "inst", "validation", "SCIENTIFIC-LIMITATIONS.md"),
    file.path(root, "inst", "validation", "CLAIM-MATRIX.md")
  )
  paths <- unique(c(
    paths[file.exists(paths)],
    system.file("validation", "SCIENTIFIC-LIMITATIONS.md", package = "bhfvar")
  ))
  paths <- paths[nzchar(paths) & file.exists(paths)]
  if (!length(paths)) stop("No public documentation files are available")
  # Collapse whitespace so that content assertions are insensitive to how the
  # documentation happens to be line-wrapped.
  gsub("\\s+", " ", paste(vapply(paths, function(path) {
    paste(readLines(path, warn = FALSE), collapse = " ")
  }, character(1)), collapse = " "))
}

# The migration vignette documents the *previous* API on purpose, so it is
# excluded from the stale-API guard below. Path resolution mirrors
# read_public_text(), including the installed-package fallback, so the guard
# still has documentation to read under `R CMD check`.
read_public_text_excluding_migration <- function() {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  paths <- c(
    file.path(root, "DESCRIPTION"),
    file.path(root, "README.md"),
    file.path(root, "NEWS.md"),
    file.path(root, "R", "bhfvar-package.R"),
    list.files(file.path(root, "vignettes"), pattern = "[.]Rmd$",
               full.names = TRUE),
    file.path(root, "inst", "validation", "SCIENTIFIC-LIMITATIONS.md"),
    file.path(root, "inst", "validation", "CLAIM-MATRIX.md")
  )
  paths <- unique(c(
    paths[file.exists(paths)],
    system.file("validation", "SCIENTIFIC-LIMITATIONS.md", package = "bhfvar")
  ))
  paths <- paths[nzchar(paths) & file.exists(paths)]
  paths <- paths[!grepl("migration", basename(paths), fixed = TRUE)]
  if (!length(paths)) stop("No public documentation files are available")
  # Collapse whitespace so that content assertions are insensitive to how the
  # documentation happens to be line-wrapped.
  gsub("\\s+", " ", paste(vapply(paths, function(path) {
    paste(readLines(path, warn = FALSE), collapse = " ")
  }, character(1)), collapse = " "))
}

test_that("public documentation gives one scoped coverage record", {
  text <- read_public_text()
  expect_match(text, "not guaranteed", ignore.case = TRUE)
  expect_match(text, "pseudo-posterior credible interval", ignore.case = TRUE)
  expect_match(text, "0.900", fixed = TRUE)
  expect_match(text, "0.871", fixed = TRUE)
  expect_match(text, "0.928", fixed = TRUE)
  expect_match(text, "0.860", fixed = TRUE)
  expect_match(text, "0.930", fixed = TRUE)
  expect_match(text, "5 percentage-point|5 percentage point",
               ignore.case = TRUE)
  expect_match(text, "reduced article-aligned balanced synthetic",
               ignore.case = TRUE)
  expect_match(text, "not propagated")
  # The article's restricted-data application must be disclaimed.
  expect_match(text, "restricted", ignore.case = TRUE)
  expect_match(text, "not reproduce", ignore.case = TRUE)
})

test_that("scientific interpretation and article locations are constrained", {
  text <- read_public_text()
  expect_match(text, "Appendix B", fixed = TRUE)
  expect_match(text, "design-effect-zero reference standardization",
               ignore.case = TRUE)
  expect_match(text, "standardization-sensitivity", ignore.case = TRUE)
  expect_match(text, "pre-centering", ignore.case = TRUE)
  expect_match(text, "boundary-dominated", ignore.case = TRUE)
  expect_match(text, "New submission", fixed = TRUE)

  banned <- c(
    "model specified in Appendix A", "model given in Appendix A",
    "implements the model in the article's Appendix A",
    "(policy estimand)", "diagnostic for informative sampling",
    "article's design family", "design used in the accompanying article",
    "coverage was still below nominal", "coverage remained below nominal",
    "Neither is a way out", "neither is a way out",
    "Gaps below roughly 1.5%", "A* is close to being undefined",
    "no errors, warnings, or notes"
  )
  for (term in banned) {
    expect_false(grepl(term, text, fixed = TRUE),
                 info = paste("banned scientific claim present:", term))
  }
})

test_that("internal governance vocabulary is absent from public documentation", {
  text <- read_public_text()
  banned <- c("Gate S4", "G4-R1", "G4-R0", "Gate E4", "Gate E5",
              "19/24", "local engineering candidate", "recovery gate")
  for (term in banned) {
    expect_false(grepl(term, text, fixed = TRUE),
                 info = paste("banned internal term present:", term))
  }
})

test_that("stale public API and visualization claims are absent", {
  text <- read_public_text_excluding_migration()
  expect_false(grepl("fit\\$diagnostics\\$ess_min", text))
  expect_false(grepl("control = list\\(", text))
  expect_false(grepl("domain_estimates\\(fit, type = \\\"marginal\\\"", text))
  expect_false(grepl("posterior predictive checks \\(requires bayesplot\\)",
                     text, ignore.case = TRUE))
  # The plotting API exists as of 0.6.0; the documentation must name it and say
  # that ggplot2 is only suggested, so a reader knows why a plot call may fail.
  expect_match(text, "bhf_plot_", fixed = TRUE)
  expect_match(text, "ggplot2", fixed = TRUE)
  expect_match(text, "suggested dependency", ignore.case = TRUE)
})

test_that("documented public exports match the namespace contract", {
  expected <- c(
    # core workflow
    "bhf_fit", "calc_eff_n", "compile_bhf_model", "domain_estimates",
    "get_stan_file_path", "log_lik", "overall_estimate",
    "prepare_bhf_data", "variance_decomposition",
    # plotting (ggplot2 is a Suggests dependency)
    "bhf_plot_variance", "bhf_plot_icc", "bhf_plot_domains",
    "bhf_plot_shrinkage", "bhf_plot_astar_sensitivity",
    # diagnostics
    "bhf_astar_sensitivity"
  )
  expect_setequal(getNamespaceExports("bhfvar"), expected)

  expect_identical(names(formals(variance_decomposition)),
                   c("fit", "prob", "print"))
  expect_identical(names(formals(domain_estimates)),
                   c("fit", "estimand", "prob", "type"))
  expect_identical(names(formals(overall_estimate)),
                   c("fit", "estimand", "prob"))
  expect_identical(names(formals(log_lik)), c("fit", "kind", "aggregate"))
})

test_that("every vignette declares an execution policy", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  vignette_dir <- file.path(root, "vignettes")
  if (!dir.exists(vignette_dir)) {
    succeed("source vignettes are validated before archive isolation")
    return(invisible(NULL))
  }
  paths <- list.files(vignette_dir, pattern = "[.]Rmd$", full.names = TRUE)
  text <- vapply(paths, function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
  }, character(1))
  expect_true(all(grepl("knitr::opts_chunk\\$set", text)))
  expect_true(any(grepl("\\{r prepare\\}", text)))
})
