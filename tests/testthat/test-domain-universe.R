test_that("supplied population shares align exactly by observed domain name", {
  result <- resolve_population_shares(
    population_shares = c(C = 5, A = 2, B = 3),
    domain_labels = c("A", "B", "C"),
    state_id = c(1L, 2L, 3L),
    raw_weights = c(1, 1, 1)
  )

  expect_s3_class(result, "bhf_population_shares")
  expect_identical(names(result$values), c("A", "B", "C"))
  expect_equal(unname(result$values), c(0.2, 0.3, 0.5))
  expect_equal(sum(result$values), 1)
  expect_identical(result$source, "external_known")
  expect_equal(result$input_sum, 10)
  expect_equal(result$normalization_factor, 0.1)
})

test_that("population share names must be present, unique, and non-blank", {
  common <- list(
    domain_labels = c("A", "B"),
    state_id = c(1L, 2L),
    raw_weights = c(1, 1)
  )

  expect_error(
    do.call(resolve_population_shares, c(
      list(population_shares = c(0.4, 0.6)), common
    )),
    "must be a named vector",
    fixed = TRUE
  )
  expect_error(
    do.call(resolve_population_shares, c(
      list(population_shares = c(A = 0.4, A = 0.6)), common
    )),
    "names must be unique",
    fixed = TRUE
  )
  expect_error(
    do.call(resolve_population_shares, c(
      list(population_shares = stats::setNames(c(0.4, 0.6), c("A", " "))),
      common
    )),
    "names must be non-missing and non-blank",
    fixed = TRUE
  )
  expect_error(
    do.call(resolve_population_shares, c(
      list(population_shares = stats::setNames(c(0.4, 0.6), c("A", NA))),
      common
    )),
    "names must be non-missing and non-blank",
    fixed = TRUE
  )
})

test_that("population share names must equal the observed domain set", {
  expect_error(
    resolve_population_shares(
      population_shares = c(A = 0.4, C = 0.6),
      domain_labels = c("A", "B"),
      state_id = c(1L, 2L),
      raw_weights = c(1, 1)
    ),
    "Missing: B. Unknown: C.",
    fixed = TRUE
  )
})

test_that("population share values must be numeric, finite, and positive", {
  args <- list(
    domain_labels = c("A", "B"),
    state_id = c(1L, 2L),
    raw_weights = c(1, 1)
  )

  expect_error(
    do.call(resolve_population_shares, c(
      list(population_shares = c(A = "0.4", B = "0.6")), args
    )),
    "must be a numeric vector",
    fixed = TRUE
  )

  invalid_values <- list(
    c(A = NA_real_, B = 1),
    c(A = Inf, B = 1),
    c(A = 0, B = 1),
    c(A = -1, B = 2)
  )
  for (shares in invalid_values) {
    expect_error(
      do.call(resolve_population_shares, c(
        list(population_shares = shares), args
      )),
      "must contain only finite, strictly positive values",
      fixed = TRUE
    )
  }
})

test_that("population share normalization requires a finite sum", {
  expect_error(
    resolve_population_shares(
      population_shares = c(
        A = .Machine$double.xmax,
        B = .Machine$double.xmax
      ),
      domain_labels = c("A", "B"),
      state_id = c(1L, 2L),
      raw_weights = c(1, 1)
    ),
    "normalization sum must be finite",
    fixed = TRUE
  )

  expect_error(
    resolve_population_shares(
      population_shares = c(
        A = .Machine$double.xmin,
        B = .Machine$double.xmax
      ),
      domain_labels = c("A", "B"),
      state_id = c(1L, 2L),
      raw_weights = c(1, 1)
    ),
    "Normalized population shares must remain finite and strictly positive.",
    fixed = TRUE
  )
})

test_that("NULL shares are estimated from retained raw weights", {
  result <- resolve_population_shares(
    population_shares = NULL,
    domain_labels = c("A", "B"),
    state_id = c(1L, 2L, 2L, 1L),
    raw_weights = c(1, 2, 4, 3)
  )

  expect_identical(names(result$values), c("A", "B"))
  expect_equal(unname(result$values), c(0.4, 0.6))
  expect_identical(result$source, "estimated_from_raw_weights")
  expect_equal(result$input_sum, 10)
})

test_that("raw-row and supplied-name permutations preserve resolved shares", {
  baseline <- resolve_population_shares(
    population_shares = c(A = 2, B = 3, C = 5),
    domain_labels = c("A", "B", "C"),
    state_id = c(1L, 2L, 3L, 1L),
    raw_weights = c(1, 2, 3, 4)
  )
  permuted_names <- resolve_population_shares(
    population_shares = c(C = 5, A = 2, B = 3),
    domain_labels = c("A", "B", "C"),
    state_id = c(1L, 2L, 3L, 1L),
    raw_weights = c(1, 2, 3, 4)
  )

  order <- c(4L, 2L, 1L, 3L)
  estimated <- resolve_population_shares(
    population_shares = NULL,
    domain_labels = c("A", "B", "C"),
    state_id = c(1L, 2L, 3L, 1L),
    raw_weights = c(1, 2, 3, 4)
  )
  estimated_permuted <- resolve_population_shares(
    population_shares = NULL,
    domain_labels = c("A", "B", "C"),
    state_id = c(1L, 2L, 3L, 1L)[order],
    raw_weights = c(1, 2, 3, 4)[order]
  )

  expect_equal(permuted_names$values, baseline$values)
  expect_equal(estimated_permuted$values, estimated$values)
})

test_that("unused factor levels do not enter the observed domain universe", {
  labels <- factor(c("B", "A"), levels = c("A", "B", "UNUSED"))
  result <- resolve_population_shares(
    population_shares = c(A = 3, B = 2),
    domain_labels = labels,
    state_id = c(1L, 2L),
    raw_weights = c(1, 1)
  )

  expect_identical(names(result$values), c("B", "A"))
  expect_equal(unname(result$values), c(0.4, 0.6))
  expect_false("UNUSED" %in% names(result$values))
})

test_that("domain IDs and raw weights define a valid observed universe", {
  expect_error(
    resolve_population_shares(
      NULL, c("A", "B"), c(1L, 1L), c(1, 2)
    ),
    "must contain only observed domains",
    fixed = TRUE
  )
  expect_error(
    resolve_population_shares(
      NULL, c("A", "B"), c(1L, 3L), c(1, 2)
    ),
    "must be between 1 and `length(domain_labels)`",
    fixed = TRUE
  )
  expect_error(
    resolve_population_shares(
      NULL, c("A", "B"), c(1L, 2L), c(1, 0)
    ),
    "raw_weights` must contain only finite, strictly positive values",
    fixed = TRUE
  )
})

test_that("prepare_bhf_data aligns supplied shares by name", {
  fixture <- make_tiny_crossed_design_fixture()

  expect_warning(
    prepared <- prepare_bhf_data(
      fixture$data,
      outcome = "outcome",
      domain = "state",
      strata = "stratum",
      psu = "psu",
      weights = "weight",
      population_shares = c(C = 2, A = 5, B = 3)
    ),
    "Sample size"
  )

  expect_equal(as.vector(prepared$stan_data$w_state_pop_share),
               c(0.5, 0.3, 0.2), tolerance = 1e-12)
  expect_identical(names(prepared$population_share_info$values),
                   c("A", "B", "C"))
  expect_identical(prepared$population_share_info$source, "external_known")
  expect_identical(prepared$input_info$population_share_source,
                   "external_known")

  expect_error(
    suppressWarnings(prepare_bhf_data(
      fixture$data,
      "outcome", "state", "stratum", "psu", "weight",
      population_shares = c(0.5, 0.3, 0.2)
    )),
    "must be a named vector"
  )
})

test_that("estimated population shares use raw weights under either scaling", {
  fixture <- make_tiny_crossed_design_fixture()

  mean_one <- suppressWarnings(prepare_bhf_data(
    fixture$data,
    "outcome", "state", "stratum", "psu", "weight",
    weight_scaling = "mean_one"
  ))
  legacy <- suppressWarnings(prepare_bhf_data(
    fixture$data,
    "outcome", "state", "stratum", "psu", "weight",
    weight_scaling = "legacy_d2"
  ))
  expected <- vapply(
    c("A", "B", "C"),
    function(label) sum(fixture$data$weight[fixture$data$state == label]),
    numeric(1)
  )
  expected <- expected / sum(expected)

  expect_equal(mean_one$population_share_info$values, expected,
               tolerance = 1e-12)
  expect_equal(legacy$population_share_info$values, expected,
               tolerance = 1e-12)
  expect_identical(mean_one$population_share_info$source,
                   "estimated_from_raw_weights")
})
