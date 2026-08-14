make_taylor_vhat_fixture <- function() {
  data <- data.frame(
    stratum = rep(c("h1", "h2"), each = 8),
    raw_psu = rep(rep(c("p1", "p2"), each = 4), 2),
    domain_label = rep(c("A", "A", "B", "B"), 4),
    y = c(
      0, 1, 0, 0, 1, 1, 1, 0,
      0, 0, 1, 1, 1, 0, 0, 1
    ),
    raw_weight = c(
      1, 2, 1, 1, 2, 1, 1, 2,
      1, 1, 2, 1, 1, 2, 1, 2
    ),
    stringsAsFactors = FALSE
  )
  data$stratum_id <- match(data$stratum, c("h1", "h2"))
  data$psu_flat_id <- match(
    paste(data$stratum, data$raw_psu, sep = "::"),
    c("h1::p1", "h1::p2", "h2::p1", "h2::p2")
  )
  data
}

direct_survey_vhat <- function(data, domain_labels) {
  direct_data <- data
  direct_data$domain_label <- factor(
    direct_data$domain_label,
    levels = domain_labels
  )
  design <- survey::svydesign(
    ids = ~psu_flat_id,
    strata = ~stratum_id,
    weights = ~raw_weight,
    data = direct_data,
    nest = TRUE
  )
  estimates <- survey::svyby(
    ~y,
    ~domain_label,
    design,
    survey::svymean,
    vartype = "se",
    keep.var = TRUE,
    covmat = TRUE,
    na.rm = FALSE,
    drop.empty.groups = FALSE
  )
  labels <- as.character(estimates$domain_label)
  stats::setNames(as.double(survey::SE(estimates))^2, labels)[domain_labels]
}

test_that("Taylor adapter matches a direct survey calculation", {
  skip_if_not_installed("survey")
  data <- make_taylor_vhat_fixture()
  domain_labels <- c("A", "B")
  direct <- direct_survey_vhat(data, domain_labels)

  result <- estimate_taylor_vhat(data, domain_labels)

  expect_true(result$enabled)
  expect_true(result$a_star_available)
  expect_equal(result$named_values, direct, tolerance = 1e-12)
  expect_equal(result$named_values, c(A = 146 / 14641, B = 146 / 14641),
               tolerance = 1e-12)
  expect_identical(result$stan_values, unname(result$named_values))
  expect_identical(result$provenance$mode, "taylor")
  expect_identical(result$provenance$source, "builtin_survey_taylor")
  expect_identical(result$provenance$weights, "raw")
  expect_identical(result$provenance$lonely_psu, "fail")
  expect_identical(result$provenance$n_strata, 2L)
  expect_identical(result$provenance$n_psus, 4L)
  expect_false(result$provenance$fpc)
  expect_false(result$provenance$uncertainty_propagated)
})

test_that("Taylor adapter uses composite PSU IDs with repeated raw labels", {
  skip_if_not_installed("survey")
  data <- make_taylor_vhat_fixture()

  expect_equal(length(unique(data$raw_psu)), 2L)
  expect_equal(length(unique(data$psu_flat_id)), 4L)
  expect_equal(
    estimate_taylor_vhat(data, c("A", "B"))$named_values,
    c(A = 146 / 14641, B = 146 / 14641),
    tolerance = 1e-12
  )
})

test_that("Taylor adapter is row-permutation equivariant", {
  skip_if_not_installed("survey")
  data <- make_taylor_vhat_fixture()
  expected <- estimate_taylor_vhat(data, c("A", "B"))$named_values
  permutation <- c(16:1)

  observed <- estimate_taylor_vhat(
    data[permutation, , drop = FALSE],
    c("A", "B")
  )$named_values

  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("Taylor adapter is invariant to a global raw-weight scale", {
  skip_if_not_installed("survey")
  data <- make_taylor_vhat_fixture()
  expected <- estimate_taylor_vhat(data, c("A", "B"))$named_values

  data$raw_weight <- data$raw_weight * 1e8
  observed <- estimate_taylor_vhat(data, c("A", "B"))$named_values

  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("Taylor mode ignores unrelated likelihood weights", {
  skip_if_not_installed("survey")
  data <- make_taylor_vhat_fixture()
  data$w_lik <- data$raw_weight / mean(data$raw_weight)
  first <- estimate_taylor_vhat(data, c("A", "B"))$named_values

  data$w_lik <- rev(data$w_lik) * 100
  second <- estimate_taylor_vhat(data, c("A", "B"))$named_values

  expect_identical(second, first)
})

test_that("Taylor adapter restores survey.lonely.psu after success", {
  skip_if_not_installed("survey")
  old <- options(survey.lonely.psu = "average")
  on.exit(options(old), add = TRUE)

  estimate_taylor_vhat(make_taylor_vhat_fixture(), c("A", "B"))

  expect_identical(getOption("survey.lonely.psu"), "average")
})

test_that("Taylor adapter restores previously unset survey options", {
  skip_if_not_installed("survey")
  option_names <- c("survey.lonely.psu", "survey.adjust.domain.lonely")
  original <- options()
  was_set <- option_names %in% names(original)
  values <- lapply(option_names, function(name) original[[name]])
  on.exit({
    for (index in seq_along(option_names)) {
      value <- if (was_set[index]) values[[index]] else NULL
      options(stats::setNames(list(value), option_names[index]))
    }
  }, add = TRUE)

  for (name in option_names) {
    options(stats::setNames(list(NULL), name))
  }
  estimate_taylor_vhat(make_taylor_vhat_fixture(), c("A", "B"))
  expect_false(any(option_names %in% names(options())))

  singleton <- make_taylor_vhat_fixture()
  singleton <- singleton[!(singleton$stratum_id == 2L &
                             singleton$psu_flat_id == 4L), , drop = FALSE]
  expect_error(
    estimate_taylor_vhat(singleton, c("A", "B")),
    class = "bhf_taylor_singleton_error"
  )
  expect_false(any(option_names %in% names(options())))
})

test_that("Taylor adapter fails closed for singleton strata and restores option", {
  skip_if_not_installed("survey")
  data <- make_taylor_vhat_fixture()
  data <- data[!(data$stratum_id == 2L & data$psu_flat_id == 4L), , drop = FALSE]
  old <- options(survey.lonely.psu = "adjust")
  on.exit(options(old), add = TRUE)

  expect_error(
    estimate_taylor_vhat(data, c("A", "B")),
    "at least two PSUs per stratum",
    class = "bhf_taylor_design_error"
  )
  expect_identical(getOption("survey.lonely.psu"), "adjust")
})

test_that("Taylor adapter rejects a PSU assigned to multiple strata", {
  skip_if_not_installed("survey")
  data <- make_taylor_vhat_fixture()
  data$psu_flat_id[data$stratum_id == 2L & data$psu_flat_id == 3L] <- 1L

  expect_error(
    estimate_taylor_vhat(data, c("A", "B")),
    "exactly one stratum",
    class = "bhf_taylor_design_error"
  )
})

test_that("Taylor adapter reports a missing survey dependency", {
  local_mocked_bindings(
    .bhf_survey_available = function() FALSE,
    .package = "bhfvar"
  )

  expect_error(
    estimate_taylor_vhat(make_taylor_vhat_fixture(), c("A", "B")),
    "requires the suggested package 'survey'",
    class = "bhf_taylor_dependency_error"
  )
})

test_that("Taylor resolver returns the adapter contract", {
  skip_if_not_installed("survey")
  result <- resolve_sampling_variances(
    deattenuation = "taylor",
    domain_labels = c("A", "B"),
    analysis_data = make_taylor_vhat_fixture()
  )

  expect_true(result$enabled)
  expect_true(result$a_star_available)
  expect_identical(result$provenance$mode, "taylor")
  expect_equal(
    result$named_values,
    c(A = 146 / 14641, B = 146 / 14641),
    tolerance = 1e-12
  )
})
