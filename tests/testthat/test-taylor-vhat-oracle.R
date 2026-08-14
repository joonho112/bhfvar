# This fixture and oracle deliberately call no bhfvar function. They provide an
# independently derived with-replacement Taylor variance for two domain Hajek
# proportions under a stratified, clustered design.
.wr_taylor_gold_fixture <- function() {
  balanced_weight <- 89 / 8
  pieces <- vector("list", 8L)
  piece <- 0L

  for (domain_label in c("A", "B")) {
    for (stratum_id in seq_len(2L)) {
      for (psu_in_stratum in seq_len(2L)) {
        piece <- piece + 1L
        base_weight <- if (stratum_id == 1L) 11 else 5
        pieces[[piece]] <- data.frame(
          y = c(if (psu_in_stratum == 1L) 1 else 0, 0, 1),
          domain_label = rep.int(domain_label, 3L),
          state_id = rep.int(match(domain_label, c("A", "B")), 3L),
          stratum_id = rep.int(stratum_id, 3L),
          psu_flat_id = rep.int(
            2L * (stratum_id - 1L) + psu_in_stratum,
            3L
          ),
          raw_weight = c(base_weight, balanced_weight, balanced_weight),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, pieces)
}

.manual_wr_hajek_vhat <- function(analysis_data, domain_labels) {
  required <- c(
    "y", "domain_label", "stratum_id", "psu_flat_id", "raw_weight"
  )
  if (!all(required %in% names(analysis_data))) {
    stop("manual oracle fixture is missing a required column", call. = FALSE)
  }

  psu_frame <- unique(analysis_data[c("stratum_id", "psu_flat_id")])
  psus_per_stratum <- table(psu_frame$stratum_id)
  if (any(psus_per_stratum < 2L)) {
    stop(
      "manual WR Taylor oracle requires at least two PSUs per stratum",
      call. = FALSE
    )
  }

  estimates <- numeric(length(domain_labels))
  numerators <- numeric(length(domain_labels))
  variances <- numeric(length(domain_labels))

  for (s in seq_along(domain_labels)) {
    in_domain <- analysis_data$domain_label == domain_labels[s]
    domain_weight <- sum(analysis_data$raw_weight[in_domain])
    estimate <- sum(
      analysis_data$raw_weight[in_domain] * analysis_data$y[in_domain]
    ) / domain_weight

    linearized_total <- analysis_data$raw_weight * in_domain *
      (analysis_data$y - estimate)
    psu_totals <- stats::aggregate(
      linearized_total,
      by = list(
        stratum_id = analysis_data$stratum_id,
        psu_flat_id = analysis_data$psu_flat_id
      ),
      FUN = sum
    )

    numerator <- sum(vapply(
      split(psu_totals$x, psu_totals$stratum_id),
      function(psu_values) {
        n_psu <- length(psu_values)
        n_psu / (n_psu - 1) * sum((psu_values - mean(psu_values))^2)
      },
      numeric(1)
    ))

    estimates[s] <- estimate
    numerators[s] <- numerator
    variances[s] <- numerator / domain_weight^2
  }

  list(
    estimate = stats::setNames(estimates, domain_labels),
    numerator = stats::setNames(numerators, domain_labels),
    variance = stats::setNames(variances, domain_labels)
  )
}

.direct_survey_vhat <- function(analysis_data) {
  design <- survey::svydesign(
    ids = ~psu_flat_id,
    strata = ~stratum_id,
    weights = ~raw_weight,
    data = analysis_data,
    nest = TRUE
  )
  estimates <- survey::svyby(
    ~y,
    ~domain_label,
    design = design,
    FUN = survey::svymean,
    keep.var = TRUE
  )
  stats::setNames(
    as.numeric(survey::SE(estimates))^2,
    as.character(estimates$domain_label)
  )
}

test_that("manual WR Taylor oracle has the frozen rational gold value", {
  fixture <- .wr_taylor_gold_fixture()
  oracle <- .manual_wr_hajek_vhat(fixture, c("A", "B"))
  expected <- c(A = 146 / 14641, B = 146 / 14641)

  expect_equal(
    stats::aggregate(raw_weight ~ domain_label, fixture, sum)$raw_weight,
    c(121, 121),
    tolerance = 0
  )
  expect_identical(oracle$estimate, c(A = 0.5, B = 0.5))
  expect_equal(oracle$numerator, c(A = 146, B = 146), tolerance = 1e-14)
  expect_equal(oracle$variance, expected, tolerance = 1e-14)
})

test_that("manual, survey, and package Taylor variances agree", {
  skip_if_not_installed("survey")
  fixture <- .wr_taylor_gold_fixture()
  expected <- c(A = 146 / 14641, B = 146 / 14641)

  manual <- .manual_wr_hajek_vhat(fixture, c("A", "B"))$variance
  direct <- .direct_survey_vhat(fixture)
  package <- estimate_taylor_vhat(
    analysis_data = fixture,
    domain_labels = c("A", "B")
  )

  expect_true(package$enabled)
  expect_true(package$a_star_available)
  expect_identical(package$provenance$mode, "taylor")
  expect_identical(package$provenance$source, "builtin_survey_taylor")
  expect_true(package$provenance$fixed_input)
  expect_false(package$provenance$uncertainty_propagated)
  expect_equal(manual, expected, tolerance = 1e-12)
  expect_equal(direct[names(expected)], expected, tolerance = 1e-12)
  expect_equal(package$named_values, expected, tolerance = 1e-12)
  expect_equal(package$named_values, direct[names(expected)], tolerance = 1e-12)
  expect_identical(package$stan_values, unname(package$named_values))
})

test_that("Taylor variance is invariant to analysis-row permutation", {
  skip_if_not_installed("survey")
  fixture <- .wr_taylor_gold_fixture()
  order <- rev(seq_len(nrow(fixture)))

  baseline <- estimate_taylor_vhat(fixture, c("A", "B"))
  permuted <- estimate_taylor_vhat(fixture[order, ], c("A", "B"))

  expect_equal(permuted$named_values, baseline$named_values, tolerance = 1e-12)
})

test_that("Taylor gold values round-trip through supplied mode by name", {
  expected <- c(A = 146 / 14641, B = 146 / 14641)
  supplied <- resolve_sampling_variances(
    deattenuation = "supplied",
    sampling_variances = expected[c("B", "A")],
    sampling_variance_method = "external_taylor",
    domain_labels = c("A", "B")
  )

  expect_identical(names(supplied$named_values), c("A", "B"))
  expect_equal(supplied$named_values, expected, tolerance = 0)
  expect_identical(supplied$stan_values, unname(expected))
})

test_that("manual and package Taylor paths reject singleton strata", {
  fixture <- .wr_taylor_gold_fixture()
  singleton <- fixture[fixture$psu_flat_id != 2L, ]
  singleton$psu_flat_id <- match(
    singleton$psu_flat_id,
    sort(unique(singleton$psu_flat_id))
  )

  expect_error(
    .manual_wr_hajek_vhat(singleton, c("A", "B")),
    "requires at least two PSUs per stratum",
    fixed = TRUE
  )
  expect_error(
    estimate_taylor_vhat(singleton, c("A", "B")),
    class = "bhf_taylor_singleton_error"
  )
})

test_that("Taylor adapter fails closed on target-label mismatch", {
  fixture <- .wr_taylor_gold_fixture()

  expect_error(
    estimate_taylor_vhat(fixture, c("A", "C")),
    class = "bhf_taylor_design_error"
  )

  malformed_labels <- list(
    c("A", "A"),
    c("A", ""),
    c("A", NA_character_)
  )
  for (labels in malformed_labels) {
    expect_error(
      estimate_taylor_vhat(fixture, labels),
      class = "bhf_vhat_alignment_error"
    )
  }
})
