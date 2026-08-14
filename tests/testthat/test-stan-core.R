# Phase 4 tests intentionally stop at Stan translation.  Native compilation
# and sampling belong to the later acceptance gate, not this static contract.

.bhf_stan_core_path <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "inst", "stan", "bhf_hybrid.stan"),
    file.path("inst", "stan", "bhf_hybrid.stan"),
    system.file("stan", "bhf_hybrid.stan", package = "bhfvar")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop("Cannot locate inst/stan/bhf_hybrid.stan", call. = FALSE)
  }
  normalizePath(existing[[1L]], mustWork = TRUE)
}

.bhf_stan_core_source <- function() {
  lines <- readLines(.bhf_stan_core_path(), warn = FALSE)
  paste(sub("//.*$", "", lines), collapse = "\n")
}

.bhf_stan_core_block <- function(name) {
  source <- .bhf_stan_core_source()
  heading <- gsub(" ", "[[:space:]]+", name, fixed = TRUE)
  hit <- regexpr(
    paste0("(?m)^[[:space:]]*", heading, "[[:space:]]*\\{"),
    source,
    perl = TRUE
  )
  if (hit[[1L]] < 0L) {
    stop("Cannot locate Stan block: ", name, call. = FALSE)
  }

  matched <- regmatches(source, hit)
  open_offset <- regexpr("\\{", matched, perl = TRUE)[[1L]] - 1L
  open_position <- hit[[1L]] + open_offset
  tail_chars <- strsplit(
    substring(source, open_position),
    split = "",
    fixed = TRUE
  )[[1L]]
  depth <- cumsum((tail_chars == "{") - (tail_chars == "}"))
  close_offset <- which(depth == 0L)[[1L]]
  substring(source, open_position + 1L, open_position + close_offset - 2L)
}

.bhf_stan_compact <- function(x) {
  trimws(gsub("[[:space:]]+", " ", x))
}

.bhf_count_matches <- function(x, pattern) {
  hits <- gregexpr(pattern, x, perl = TRUE)[[1L]]
  if (identical(hits, -1L)) 0L else length(hits)
}

.bhf_prepare_core_fixture <- function() {
  fixture <- make_tiny_crossed_design_fixture()
  prepared <- suppressWarnings(prepare_bhf_data(
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
  list(fixture = fixture, prepared = prepared)
}

test_that("Stan data symbols and prepared-data fields share one versioned contract", {
  data_block <- .bhf_stan_compact(.bhf_stan_core_block("data"))
  x <- .bhf_prepare_core_fixture()
  stan_data <- x$prepared$stan_data

  expected_fields <- c(
    "data_schema_version", "N", "S", "H", "J", "y", "state_id",
    "stratum_id", "J_h", "psu_start", "psu_in_stratum_id",
    "psu_flat_id", "w_lik", "w_state_pop_share", "prior_alpha_mean",
    "prior_alpha_sd", "sigma_state_prior_code", "use_deattenuation",
    "vhat_state"
  )
  expect_setequal(names(stan_data), expected_fields)
  expect_silent(validate_stan_data(stan_data))
  expect_identical(stan_data$data_schema_version, 1L)
  expect_identical(
    c(N = stan_data$N, S = stan_data$S, H = stan_data$H, J = stan_data$J),
    x$fixture$truth$dimensions
  )
  expect_identical(length(stan_data$state_id), stan_data$N)
  expect_identical(length(stan_data$stratum_id), stan_data$N)
  expect_identical(length(stan_data$psu_flat_id), stan_data$N)
  expect_identical(length(stan_data$J_h), stan_data$H)
  expect_identical(length(stan_data$w_state_pop_share), stan_data$S)

  declarations <- c(
    "int<lower=1> data_schema_version;",
    "int<lower=1> N;", "int<lower=1> S;", "int<lower=1> H;",
    "int<lower=1> J;",
    "array[N] int<lower=0, upper=1> y;",
    "array[N] int<lower=1, upper=S> state_id;",
    "array[N] int<lower=1, upper=H> stratum_id;",
    "array[H] int<lower=1> J_h;",
    "array[H] int<lower=1> psu_start;",
    "array[N] int<lower=1> psu_in_stratum_id;",
    "array[N] int<lower=1, upper=J> psu_flat_id;",
    "vector<lower=0>[N] w_lik;",
    "vector<lower=0, upper=1>[S] w_state_pop_share;",
    "real prior_alpha_mean;", "real<lower=0> prior_alpha_sd;",
    "int<lower=1, upper=4> sigma_state_prior_code;",
    "int<lower=0, upper=1> use_deattenuation;",
    "vector<lower=0>[S] vhat_state;"
  )
  for (declaration in declarations) {
    pattern <- gsub("([][{}()+*.^$|\\\\?])", "\\\\\\1", declaration)
    pattern <- gsub(" ", "[[:space:]]*", pattern, fixed = TRUE)
    expect_match(data_block, pattern, perl = TRUE, info = declaration)
  }
})

test_that("transformed data rejects malformed indexing and transport inputs", {
  transformed_data <- .bhf_stan_compact(
    .bhf_stan_core_block("transformed data")
  )
  all_source <- .bhf_stan_compact(.bhf_stan_core_source())

  required_guards <- c(
    "data_schema_version != 1",
    "min\\(w_lik\\) <= 0",
    "abs\\(sum\\(w_lik\\) - N\\)",
    "min\\(w_state_pop_share\\) <= 0",
    "abs\\(sum\\(w_state_pop_share\\) - 1\\.0\\)",
    "prior_alpha_sd <= 0",
    "use_deattenuation == 0 && max\\(vhat_state\\) > 0",
    "sum\\(J_h\\) != J",
    "psu_start\\[1\\] != 1",
    "psu_start\\[h\\] != psu_start\\[h - 1\\] \\+ J_h\\[h - 1\\]",
    "psu_in_stratum_id\\[i\\] > J_h\\[h\\]",
    paste0(
      "psu_flat_id\\[i\\] != psu_start\\[h\\] \\+ ",
      "psu_in_stratum_id\\[i\\] - 1"
    )
  )
  for (guard in required_guards) {
    expect_match(transformed_data, guard, perl = TRUE, info = guard)
  }

  expect_false(grepl("\\bw_norm\\b|\\bw_sum\\b", all_source, perl = TRUE))
  expect_false(grepl(
    "w_lik\\s*\\[[^]]+\\]\\s*/\\s*sum\\(w_lik\\)",
    all_source,
    perl = TRUE
  ))
})

test_that("three non-centered random-effect families have fixed dimensions", {
  parameters <- .bhf_stan_compact(.bhf_stan_core_block("parameters"))
  transformed <- .bhf_stan_compact(
    .bhf_stan_core_block("transformed parameters")
  )

  for (pattern in c(
    "vector\\[S\\] z_state;", "vector\\[H\\] z_stratum;",
    "vector\\[J\\] z_psu;", "real<lower=0> sigma_state;",
    "real<lower=0> sigma_stratum;", "real<lower=0> sigma_psu;"
  )) {
    expect_match(parameters, pattern, perl = TRUE, info = pattern)
  }
  expect_false(grepl("\\bu_(state|stratum|psu)\\b", parameters, perl = TRUE))
  for (pattern in c(
    "vector\\[S\\] u_state;", "vector\\[H\\] u_stratum;",
    "vector\\[J\\] u_psu;", "vector\\[N\\] eta;"
  )) {
    expect_match(transformed, pattern, perl = TRUE, info = pattern)
  }
})

test_that("state, stratum, and within-stratum PSU effects are centered", {
  transformed <- .bhf_stan_compact(
    .bhf_stan_core_block("transformed parameters")
  )
  generated <- .bhf_stan_compact(
    .bhf_stan_core_block("generated quantities")
  )

  centering_patterns <- c(
    "state_weighted_mean = dot_product\\(w_state_pop_share, z_state\\)",
    paste0(
      "u_state = sigma_state \\* \\(z_state - ",
      "rep_vector\\(state_weighted_mean, S\\)\\)"
    ),
    paste0(
      "u_stratum = sigma_stratum \\* \\(z_stratum - ",
      "rep_vector\\(mean\\(z_stratum\\), H\\)\\)"
    ),
    "for \\(h in 1:H\\)",
    "for \\(j_local in 1:J_h\\[h\\]\\)",
    "j_flat = psu_start\\[h\\] \\+ j_local - 1",
    "block_mean /= J_h\\[h\\]",
    paste0(
      "u_psu\\[j_flat\\] = sigma_psu \\* ",
      "\\(z_psu\\[j_flat\\] - block_mean\\)"
    )
  )
  for (pattern in centering_patterns) {
    expect_match(transformed, pattern, perl = TRUE, info = pattern)
  }
  expect_match(
    generated,
    "dot_product\\(w_state_pop_share, u_state\\)",
    perl = TRUE
  )
  expect_match(generated, "mean\\(u_stratum\\)", perl = TRUE)
  expect_match(
    generated,
    "psu_centering_residual\\[h\\] = block_sum / J_h\\[h\\]",
    perl = TRUE
  )

  # J_h has lower bound one, and subtracting the one-element block mean makes
  # a singleton PSU effect exactly zero without an exception or special case.
  expect_equal(0.8 * (9 - mean(9)), 0)
})

test_that("the full crossed/nested predictor matches an independent hand calculation", {
  transformed <- .bhf_stan_compact(
    .bhf_stan_core_block("transformed parameters")
  )
  expect_match(
    transformed,
    paste0(
      "eta\\[i\\] = alpha \\+ u_state\\[state_id\\[i\\]\\] \\+ ",
      "u_stratum\\[stratum_id\\[i\\]\\] \\+ ",
      "u_psu\\[psu_flat_id\\[i\\]\\]"
    ),
    perl = TRUE
  )

  x <- .bhf_prepare_core_fixture()
  stan_data <- x$prepared$stan_data
  alpha <- -0.35
  u_state <- 0.7 * (c(1.1, -0.2, 0.7) - 0.63)
  u_stratum <- 0.4 * (c(-0.4, 1.0) - 0.3)
  u_psu <- 0.25 * c(1.5, -1.5, 0, 0)
  eta <- alpha + u_state[stan_data$state_id] +
    u_stratum[stan_data$stratum_id] + u_psu[stan_data$psu_flat_id]

  expect_equal(sum(stan_data$w_state_pop_share * u_state), 0,
               tolerance = 1e-15)
  expect_equal(mean(u_stratum), 0, tolerance = 1e-15)
  expect_equal(c(mean(u_psu[1:2]), mean(u_psu[3:4])), c(0, 0),
               tolerance = 1e-15)
  expect_equal(
    eta,
    c(
      0.074, -0.836, -0.206, -0.676, -1.586, -0.956,
      0.259, -0.651, -0.021, 0.259, -0.651, -0.021
    ),
    tolerance = 1e-15
  )
})

test_that("baseline priors and the state-only sensitivity selector are closed", {
  model <- .bhf_stan_compact(.bhf_stan_core_block("model"))

  required_priors <- c(
    "alpha ~ normal\\(prior_alpha_mean, prior_alpha_sd\\)",
    "z_state ~ std_normal\\(\\)",
    "z_stratum ~ std_normal\\(\\)",
    "z_psu ~ std_normal\\(\\)",
    "sigma_state_prior_code == 1",
    "sigma_state ~ student_t\\(3, 0, 2\\.5\\)",
    "sigma_state_prior_code == 2",
    "sigma_state ~ normal\\(0, 1\\)",
    "sigma_state_prior_code == 3",
    "sigma_state ~ cauchy\\(0, 2\\.5\\)",
    "sigma_state ~ student_t\\(3, 0, 5\\)",
    "sigma_stratum ~ student_t\\(3, 0, 2\\.5\\)",
    "sigma_psu ~ student_t\\(3, 0, 2\\.5\\)"
  )
  for (prior in required_priors) {
    expect_match(model, prior, perl = TRUE, info = prior)
  }
  expect_identical(.bhf_count_matches(model, "sigma_stratum\\s*~"), 1L)
  expect_identical(.bhf_count_matches(model, "sigma_psu\\s*~"), 1L)

  x <- .bhf_prepare_core_fixture()
  prior <- x$prepared$prior_info
  expect_equal(x$prepared$stan_data$prior_alpha_mean,
               stats::qlogis(sum(x$fixture$data$weight *
                                   x$fixture$data$outcome) /
                               sum(x$fixture$data$weight)),
               tolerance = 1e-15)
  expect_identical(x$prepared$stan_data$prior_alpha_sd, 0.5)
  expect_identical(x$prepared$stan_data$sigma_state_prior_code, 1L)
  expect_identical(prior$random_effect_sd$family, "half_student_t")
  expect_identical(prior$random_effect_sd$applies_to,
                   c("state", "stratum", "psu"))
  expect_true(prior$sigma_state$baseline)
})

test_that("pseudo-likelihood consumes canonical weights exactly once", {
  transformed_data <- .bhf_stan_compact(
    .bhf_stan_core_block("transformed data")
  )
  model <- .bhf_stan_compact(.bhf_stan_core_block("model"))
  source <- .bhf_stan_compact(.bhf_stan_core_source())

  weighted_target <- paste0(
    "target \\+= w_lik\\[i\\] \\* ",
    "bernoulli_logit_lpmf\\(y\\[i\\] \\| eta\\[i\\]\\)"
  )
  expect_match(model, weighted_target, perl = TRUE)
  expect_identical(.bhf_count_matches(model, "target\\s*\\+="), 1L)
  expect_match(transformed_data, "abs\\(sum\\(w_lik\\) - N\\)",
               perl = TRUE)
  expect_false(grepl("target \\+= bernoulli_logit_lpmf", model, perl = TRUE))
  expect_false(grepl("\\bw_norm\\b|\\bw_sum\\b", source, perl = TRUE))

  x <- .bhf_prepare_core_fixture()
  stan_data <- x$prepared$stan_data
  eta <- rep(c(-0.8, 0.1, 0.9), length.out = stan_data$N)
  log_prob <- stats::dbinom(stan_data$y, size = 1,
                            prob = stats::plogis(eta), log = TRUE)
  pseudo_target <- sum(stan_data$w_lik * log_prob)
  unweighted_target <- sum(log_prob)

  expect_equal(sum(stan_data$w_lik), stan_data$N, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(pseudo_target, unweighted_target,
                                tolerance = 1e-12)))
  expect_equal(pseudo_target, sum(x$prepared$analysis_data$w_lik * log_prob),
               tolerance = 1e-15)
})

test_that("legacy estimand shortcuts are absent from the Phase 4 core", {
  source <- .bhf_stan_compact(.bhf_stan_core_source())
  forbidden <- c(
    "\\bc_factor\\b", "0\\.346", "fmax\\(0\\.001",
    "\\bw_norm\\b", "\\bw_sum\\b", "\\breliability_state\\b",
    "\\bicc_prob\\b", "\\bicc_deatten\\b"
  )
  for (pattern in forbidden) {
    expect_false(grepl(pattern, source, perl = TRUE), info = pattern)
  }
})

test_that("rstan translates the Phase 4 Stan core when available", {
  skip_if_not_installed("rstan")
  translated <- NULL
  expect_no_error(
    translated <- rstan::stanc(
      file = .bhf_stan_core_path(),
      model_name = "bhf_hybrid_phase4_static_test",
      allow_undefined = FALSE,
      verbose = FALSE
    )
  )
  expect_true(isTRUE(translated$status))
  expect_true(is.character(translated$cppcode))
  expect_true(nzchar(translated$cppcode))
})
