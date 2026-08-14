make_prior_sensitivity_object <- function(variant = "half_t3_2.5") {
  fixture <- make_tiny_crossed_design_fixture()
  suppressWarnings(prepare_bhf_data(
    fixture$data,
    "outcome", "state", "stratum", "psu", "weight",
    population_shares = fixture$known$population_shares,
    deattenuation = "supplied",
    sampling_variances = fixture$known$vhat,
    sampling_variance_method = "external_taylor",
    prior_alpha_mean = -1.5,
    prior_alpha_sd = 0.5,
    sigma_state_prior = variant
  ))
}

test_that("sigma_state prior catalog has exactly four article variants", {
  catalog <- sigma_state_prior_catalog()

  expect_identical(
    names(catalog),
    c("half_t3_2.5", "half_normal_1", "half_cauchy_2.5", "half_t3_5")
  )
  expect_identical(unname(vapply(catalog, `[[`, integer(1), "code")), 1:4)
  expect_identical(
    unname(vapply(catalog, `[[`, character(1), "family")),
    c("half_student_t", "half_normal", "half_cauchy", "half_student_t")
  )
  expect_equal(unname(vapply(catalog, `[[`, numeric(1), "scale")),
               c(2.5, 1, 2.5, 5))
})

test_that("baseline is the default and variants transport exact codes", {
  default <- make_prior_sensitivity_object()
  variants <- names(sigma_state_prior_catalog())
  objects <- lapply(variants, make_prior_sensitivity_object)

  expect_identical(default$prior_info$sigma_state$variant, "half_t3_2.5")
  expect_identical(default$prior_info$sigma_state$code, 1L)
  expect_identical(default$prior_info$sigma_state$source, "default_baseline")
  expect_true(default$prior_info$sigma_state$baseline)

  expect_identical(
    vapply(objects, function(x) x$stan_data$sigma_state_prior_code, integer(1)),
    1:4
  )
  expect_identical(
    vapply(objects, function(x) x$prior_info$sigma_state$variant, character(1)),
    variants
  )
  expect_identical(
    vapply(objects, function(x) x$input_info$sigma_state_prior, character(1)),
    variants
  )
  expect_true(all(vapply(
    objects,
    function(x) x$prior_info$sigma_state$other_random_effect_priors_fixed,
    logical(1)
  )))
  for (object in objects) {
    expect_invisible(validate_bhf_data_contract(object))
  }
})

test_that("prior sensitivity changes only the state selector", {
  variants <- names(sigma_state_prior_catalog())
  objects <- lapply(variants, make_prior_sensitivity_object)
  baseline <- objects[[1L]]

  for (object in objects[-1L]) {
    expect_identical(object$prior_info$alpha, baseline$prior_info$alpha)
    expect_identical(
      object$prior_info$random_effect_sd,
      baseline$prior_info$random_effect_sd
    )

    stan_baseline <- baseline$stan_data
    stan_variant <- object$stan_data
    stan_baseline$sigma_state_prior_code <- NULL
    stan_variant$sigma_state_prior_code <- NULL
    expect_identical(stan_variant, stan_baseline)

    expect_identical(
      object$prior_info$sigma_state$varied_component,
      "state"
    )
    expect_false(object$prior_info$sigma_state$baseline)
    expect_identical(
      object$prior_info$sigma_state$source,
      "article_sensitivity_selector"
    )
  }
})

test_that("sigma_state prior selector rejects non-catalog inputs", {
  fixture <- make_tiny_crossed_design_fixture()
  base_args <- list(
    data = fixture$data,
    outcome = "outcome",
    domain = "state",
    strata = "stratum",
    psu = "psu",
    weights = "weight",
    population_shares = fixture$known$population_shares,
    deattenuation = "none",
    prior_alpha_mean = -1.5
  )

  invalid <- list(
    "half_t3",
    c("half_t3_2.5", "half_normal_1"),
    NA_character_,
    1L
  )
  for (value in invalid) {
    args <- base_args
    args$sigma_state_prior <- value
    expect_error(
      suppressWarnings(do.call(prepare_bhf_data, args)),
      "must be exactly one of"
    )
  }
})

test_that("local and whole-object gates reject corrupted prior codes", {
  canonical <- make_prior_sensitivity_object("half_cauchy_2.5")

  broken_stan <- canonical$stan_data
  broken_stan$sigma_state_prior_code <- 5L
  expect_error(validate_stan_data(broken_stan), "integer in \\[1, 4\\]")

  broken <- canonical
  broken$stan_data$sigma_state_prior_code <- 2L
  error <- tryCatch(validate_bhf_data_contract(broken), error = identity)
  expect_s3_class(error, "bhf_data_contract_error")
  expect_identical(error$invariant, "prior.state_selector")

  broken <- canonical
  broken$prior_info$sigma_state$scale <- 99
  broken$provenance$prior <- broken$prior_info
  error <- tryCatch(validate_bhf_data_contract(broken), error = identity)
  expect_s3_class(error, "bhf_data_contract_error")
  expect_identical(error$invariant, "prior.state_selector")
})
