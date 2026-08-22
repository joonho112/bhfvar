# Runtime contract tests for plots and the A* sensitivity diagnostic.

skip_if_no_ggplot2 <- function() {
  testthat::skip_if_not_installed("ggplot2")
}

test_that("plot entry points reject objects they cannot use", {
  skip_if_no_ggplot2()
  for (f in list(bhf_plot_variance, bhf_plot_icc)) {
    expect_error(f(42), class = "bhf_argument_error")
    expect_error(f("not a fit"), class = "bhf_argument_error")
  }
  expect_error(bhf_plot_domains(42), class = "bhf_fit_class_error")
  expect_error(bhf_plot_shrinkage(42), class = "bhf_fit_class_error")
  expect_error(bhf_plot_astar_sensitivity(data.frame(x = 1)),
               class = "bhf_argument_error")
})

test_that("missing ggplot2 produces an actionable classed error", {
  local_mocked_bindings(
    .bhf_ggplot2_available = function() FALSE,
    .package = "bhfvar"
  )
  error <- expect_error(
    bhfvar:::.bhf_need_ggplot2("bhf_plot_domains"),
    class = "bhf_missing_suggest"
  )
  expect_match(conditionMessage(error), "install.packages")
  expect_match(conditionMessage(error), "tidy extractor output")
})

test_that("fit plots forward interval probability and omit unavailable A-star", {
  skip_if_no_ggplot2()
  local_result_draw_backend()

  enabled <- make_result_contract_fit(TRUE)
  variance_plot <- bhf_plot_variance(enabled, prob = 0.90)
  icc_plot <- bhf_plot_icc(enabled, prob = 0.80)
  expect_s3_class(variance_plot, "ggplot")
  expect_s3_class(icc_plot, "ggplot")
  expect_true(all(variance_plot$data$prob == 0.90))
  expect_true(all(icc_plot$data$prob == 0.80))
  expect_match(variance_plot$labels$caption, "90%")
  expect_match(icc_plot$labels$caption, "80%")

  disabled <- make_result_contract_fit(FALSE)
  no_astar <- bhf_plot_variance(disabled, prob = 0.90)
  expect_false(any(as.character(no_astar$data$estimand) == "A*"))
  expect_true(all(is.finite(no_astar$data$mean)))
})

test_that("domain plots expose only supported estimands and keep exact caps", {
  skip_if_no_ggplot2()
  local_result_draw_backend()
  fit <- make_result_contract_fit(TRUE)

  capped <- bhf_plot_domains(fit, estimand = "B", prob = 0.90,
                             n_domains = 3)
  expect_s3_class(capped, "ggplot")
  expect_equal(nrow(capped$data), 3)
  expect_true(all(capped$data$estimand == "B"))
  expect_match(capped$labels$caption, "90%")
  expect_equal(nrow(bhf_plot_domains(fit, n_domains = 2)$data), 2)
  expect_equal(nrow(bhf_plot_domains(fit, n_domains = 20)$data), 3)

  expect_error(bhf_plot_domains(fit, estimand = "A_star"),
               class = "bhf_argument_error")
  for (bad in list(1, 2.5, Inf, NA_real_, c(2, 3), "3")) {
    expect_error(bhf_plot_domains(fit, n_domains = bad),
                 class = "bhf_argument_error")
  }
})

test_that("shrinkage plot uses posterior intervals without causal wording", {
  skip_if_no_ggplot2()
  local_result_draw_backend()
  fit <- make_result_contract_fit(TRUE)

  plot <- bhf_plot_shrinkage(fit, estimand = "A", prob = 0.90)
  built <- ggplot2::ggplot_build(plot)
  expect_s3_class(plot, "ggplot")
  expect_true(any(vapply(
    built$data,
    function(layer) all(c("ymin", "ymax") %in% names(layer)),
    logical(1)
  )))
  expect_match(plot$labels$caption, "90%")
  expect_match(plot$labels$caption, "partial pooling")
  expect_false(grepl("being pulled toward", plot$labels$caption, fixed = TRUE))
  expect_error(bhf_plot_shrinkage(fit, estimand = "A_star"),
               class = "bhf_argument_error")
})

test_that("A-star sensitivity calls the result contract and estimand formula", {
  local_result_draw_backend()
  fit <- make_result_contract_fit(TRUE)
  result <- bhf_astar_sensitivity(fit, scale = c(0, 1, 2), prob = 0.90)

  expect_s3_class(result, "bhf_astar_sensitivity")
  base <- sum(c(0.5, 0.3, 0.2) * c(0.004, 0.005, 0.0075))
  expected <- vapply(c(0, 1, 2), function(s) {
    mean(pmax(0, fit$.draws$var_between_A - s * base))
  }, numeric(1))
  expect_equal(result$between_mean, expected)
  expect_equal(attr(result, "base_correction"), base)
  expect_identical(attr(result, "vhat_source"), "supplied")
  expect_true(all(diff(result$between_mean) <= 0))
  expect_true(all(diff(result$share_at_boundary) >= 0))
  expect_true(all(result$prob == 0.90))

  decomposition <- variance_decomposition(fit, prob = 0.90)
  scale_one <- result[result$scale == 1, ]
  expect_equal(
    scale_one$between_mean,
    subset(decomposition$A_star$summary, component == "between")$mean
  )
  expect_equal(
    scale_one$proportion_mean,
    subset(decomposition$A_star$summary, component == "proportion")$mean
  )

  disabled <- make_result_contract_fit(FALSE)
  expect_error(bhf_astar_sensitivity(disabled), class = "bhf_result_error")
  expect_error(bhf_astar_sensitivity(42), class = "bhf_fit_class_error")
  expect_error(bhf_astar_sensitivity(fit, scale = c(1, Inf)),
               class = "bhf_argument_error")

  missing <- fit
  missing$.draws$var_within_A <- NULL
  expect_error(bhf_astar_sensitivity(missing),
               class = "bhf_draw_contract_error")
})

test_that("plot helpers remain exported with stable signatures", {
  exported <- getNamespaceExports("bhfvar")
  for (f in c("bhf_plot_variance", "bhf_plot_icc", "bhf_plot_domains",
              "bhf_plot_shrinkage", "bhf_astar_sensitivity",
              "bhf_plot_astar_sensitivity")) {
    expect_true(f %in% exported, info = f)
  }
  expect_identical(names(formals(bhf_astar_sensitivity)),
                   c("fit", "scale", "prob"))
  expect_identical(names(formals(bhf_plot_domains)),
                   c("fit", "estimand", "prob", "n_domains"))
})
