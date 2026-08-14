.bhf_run_stan_parity <- function() {
  tolower(Sys.getenv("BHFVAR_RUN_STAN_PARITY", "false")) %in%
    c("1", "true", "yes")
}

.bhf_stan_parity_path <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "inst", "stan", "bhf_hybrid.stan"),
    file.path("inst", "stan", "bhf_hybrid.stan"),
    system.file("stan", "bhf_hybrid.stan", package = "bhfvar")
  )
  existing <- unique(candidates[nzchar(candidates) & file.exists(candidates)])
  if (!length(existing)) stop("Cannot locate Stan source", call. = FALSE)
  normalizePath(existing[[1L]], mustWork = TRUE)
}

.bhf_expect_close <- function(actual, expected, label,
                              tolerance = 1e-9) {
  normalize_scalar_draws <- function(x) {
    if (!is.null(dim(x)) && length(dim(x)) == 1L) as.vector(x) else x
  }
  expect_equal(
    unname(normalize_scalar_draws(actual)),
    unname(normalize_scalar_draws(expected)), tolerance = tolerance,
    info = label
  )
}

test_that("native parity tier is pinned to the frozen oracle artifacts", {
  files <- testthat::test_path(
    "..", c(
      "oracle/estimand_oracle.R", "oracle/smoke_oracle.R",
      "testthat/test-estimand-oracle.R", "oracle/manifest.txt"
    )
  )
  expected <- c(
    "03ff2fb87e1ca1124d671607774394b31da09a9254b79da8fd4cedf035897d72",
    "80c81b621ec0e780fa9eb7d4a17bee2eac71d8f5d7b4dc564e07e5a70d652be8",
    "ab8b54f6c20c819053aa81a8b815938051c93b0115b36a1ba376a002454d3dd9",
    "885a44a0e2fa92ca028b15301c4c17202031afbcd55504c206bd9058157dbb8a"
  )
  observed <- unname(tools::sha256sum(files))
  expect_identical(observed, expected)
})

test_that("Stan generated quantities match the frozen independent oracle", {
  if (!.bhf_run_stan_parity()) {
    succeed("Set BHFVAR_RUN_STAN_PARITY=true for the native parity tier")
  } else {
    skip_if_not_installed("rstan")

    oracle_env <- new.env(parent = baseenv())
    sys.source(
      testthat::test_path("..", "oracle", "estimand_oracle.R"),
      envir = oracle_env
    )
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
      sampling_variance_method = "external_taylor",
      prior_alpha_mean = -1.5
    ))

    model <- rstan::stan_model(
      file = .bhf_stan_parity_path(), auto_write = FALSE, verbose = FALSE
    )
    fit <- suppressWarnings(rstan::sampling(
      model,
      data = prepared$stan_data,
      chains = 1,
      iter = 180,
      warmup = 90,
      seed = 520514,
      refresh = 0,
      control = list(adapt_delta = 0.95)
    ))
    output_names <- c(
      "alpha", "u_state", "u_stratum", "u_psu", "eta",
      "sigma_state", "sigma_stratum", "sigma_psu",
      "var_state_latent", "var_stratum_latent", "var_psu_latent",
      "var_level1_latent", "var_total_latent", "icc_state_latent",
      "icc_stratum_latent", "icc_psu_latent",
      "p_state_A", "p_bar_A", "var_between_A", "var_within_A",
      "var_total_A", "prop_between_A", "A_proportion_defined",
      "mean_vhat_state", "p_bar_A_star", "var_between_A_star",
      "var_within_A_star", "var_total_A_star", "prop_between_A_star",
      "A_star_at_boundary", "A_star_truncated",
      "A_star_proportion_defined", "p_individual_B", "state_weight_B",
      "p_state_B", "p_bar_B", "within_binomial_state_B",
      "within_mixture_state_B", "var_within_binomial_B",
      "var_within_mixture_B", "var_within_B", "var_between_B",
      "var_total_B", "prop_between_B", "B_proportion_defined",
      "gap_B_minus_A_mean", "gap_B_minus_A_between",
      "gap_B_minus_A_within", "gap_B_minus_A_total",
      "gap_B_minus_A_proportion", "gap_A_minus_A_star_between",
      "gap_A_minus_A_star_total", "gap_A_minus_A_star_proportion",
      "sd_state_logit_unweighted", "sd_state_logit_weighted",
      "sd_state_probability_A_unweighted",
      "sd_state_probability_A_weighted", "log_lik_pseudo",
      "log_lik_raw", "y_rep"
    )
    draws <- rstan::extract(fit, pars = output_names, permuted = TRUE)
    expect_setequal(names(draws), output_names)

    draw_id <- as.character(seq_along(draws$alpha))
    domain_labels <- as.character(prepared$mapping$domain$label)
    stratum_labels <- as.character(prepared$mapping$stratum$label)
    psu_labels <- paste(
      prepared$mapping$psu$stratum_label,
      prepared$mapping$psu$label,
      sep = "::"
    )
    dimnames(draws$u_state) <- list(draw_id, domain_labels)
    dimnames(draws$u_stratum) <- list(draw_id, stratum_labels)
    dimnames(draws$u_psu) <- list(draw_id, psu_labels)

    design <- list(
      domain_labels = domain_labels,
      stratum_labels = stratum_labels,
      psu_labels = psu_labels,
      row_id = as.character(prepared$analysis_data$original_row),
      state_id = as.integer(prepared$stan_data$state_id),
      stratum_id = as.integer(prepared$stan_data$stratum_id),
      psu_flat_id = as.integer(prepared$stan_data$psu_flat_id),
      w_lik = as.numeric(prepared$stan_data$w_lik),
      population_share = stats::setNames(
        as.numeric(prepared$stan_data$w_state_pop_share), domain_labels
      )
    )
    oracle_draws <- list(
      draw_id = draw_id,
      alpha = draws$alpha,
      u_state = draws$u_state,
      u_stratum = draws$u_stratum,
      u_psu = draws$u_psu,
      sigma_state = draws$sigma_state,
      sigma_stratum = draws$sigma_stratum,
      sigma_psu = draws$sigma_psu
    )
    vhat <- stats::setNames(
      as.numeric(prepared$stan_data$vhat_state), domain_labels
    )
    oracle <- oracle_env$bhf_reference_oracle(
      design, oracle_draws, vhat = vhat, return_individual = TRUE
    )

    latent_actual <- cbind(
      var_state = draws$var_state_latent,
      var_stratum = draws$var_stratum_latent,
      var_psu = draws$var_psu_latent,
      var_logistic = draws$var_level1_latent,
      total = draws$var_total_latent,
      icc_state = draws$icc_state_latent,
      icc_stratum = draws$icc_stratum_latent,
      icc_psu = draws$icc_psu_latent
    )
    .bhf_expect_close(latent_actual, oracle$latent[, colnames(latent_actual)],
                      "latent diagnostics")

    a_actual <- cbind(
      mean = draws$p_bar_A,
      between = draws$var_between_A,
      within = draws$var_within_A,
      total = draws$var_total_A,
      proportion = draws$prop_between_A
    )
    .bhf_expect_close(draws$p_state_A, oracle$A$p_state, "A probabilities")
    .bhf_expect_close(a_actual, oracle$A$summary, "A decomposition")
    expect_identical(
      as.logical(draws$A_proportion_defined),
      as.logical(oracle$A$proportion_defined)
    )

    a_star_actual <- cbind(
      mean = draws$p_bar_A_star,
      between = draws$var_between_A_star,
      within = draws$var_within_A_star,
      total = draws$var_total_A_star,
      proportion = draws$prop_between_A_star
    )
    .bhf_expect_close(a_star_actual, oracle$A_star$summary,
                      "A-star decomposition")
    .bhf_expect_close(draws$mean_vhat_state, oracle$A_star$correction,
                      "A-star correction")
    expect_identical(as.logical(draws$A_star_at_boundary),
                     as.logical(oracle$A_star$at_boundary))
    expect_identical(as.logical(draws$A_star_truncated),
                     as.logical(oracle$A_star$truncated))
    expect_identical(as.logical(draws$A_star_proportion_defined),
                     as.logical(oracle$A_star$proportion_defined))

    b_actual <- cbind(
      mean = draws$p_bar_B,
      between = draws$var_between_B,
      within = draws$var_within_B,
      total = draws$var_total_B,
      proportion = draws$prop_between_B
    )
    .bhf_expect_close(draws$p_individual_B, oracle$B$p_individual,
                      "B individual probabilities")
    expected_state_weight <- vapply(
      seq_len(prepared$stan_data$S),
      function(s) sum(prepared$stan_data$w_lik[
        prepared$stan_data$state_id == s
      ]),
      numeric(1)
    )
    .bhf_expect_close(
      draws$state_weight_B,
      matrix(expected_state_weight, nrow(draws$state_weight_B),
             length(expected_state_weight), byrow = TRUE),
      "B state likelihood-weight denominators"
    )
    .bhf_expect_close(draws$p_state_B, oracle$B$p_state,
                      "B state probabilities")
    .bhf_expect_close(draws$within_binomial_state_B,
                      oracle$B$within_bernoulli_state,
                      "B state binomial terms")
    .bhf_expect_close(draws$within_mixture_state_B,
                      oracle$B$within_mixture_state,
                      "B state mixture terms")
    .bhf_expect_close(draws$var_within_binomial_B,
                      oracle$B$within_bernoulli,
                      "B aggregate binomial term")
    .bhf_expect_close(draws$var_within_mixture_B,
                      oracle$B$within_mixture,
                      "B aggregate mixture term")
    .bhf_expect_close(b_actual, oracle$B$summary, "B decomposition")
    expect_identical(as.logical(draws$B_proportion_defined),
                     as.logical(oracle$B$proportion_defined))

    .bhf_expect_close(draws$gap_B_minus_A_mean,
                      oracle$gaps$B_minus_A[, "mean"], "B-A mean gap")
    .bhf_expect_close(draws$gap_B_minus_A_between,
                      oracle$gaps$B_minus_A[, "between"], "B-A between gap")
    .bhf_expect_close(draws$gap_B_minus_A_within,
                      oracle$gaps$B_minus_A[, "within"], "B-A within gap")
    .bhf_expect_close(draws$gap_B_minus_A_total,
                      oracle$gaps$B_minus_A[, "total"], "B-A total gap")
    .bhf_expect_close(draws$gap_B_minus_A_proportion,
                      oracle$gaps$B_minus_A[, "proportion"],
                      "B-A proportion gap")
    .bhf_expect_close(draws$gap_A_minus_A_star_between,
                      oracle$gaps$A_minus_A_star[, "between"],
                      "A-A-star between gap")
    .bhf_expect_close(draws$gap_A_minus_A_star_total,
                      oracle$gaps$A_minus_A_star[, "total"],
                      "A-A-star total gap")
    .bhf_expect_close(draws$gap_A_minus_A_star_proportion,
                      oracle$gaps$A_minus_A_star[, "proportion"],
                      "A-A-star proportion gap")

    pi_state <- as.numeric(prepared$stan_data$w_state_pop_share)
    weighted_u_mean <- as.vector(draws$u_state %*% pi_state)
    sd_logit_w <- sqrt(rowSums(
      (draws$u_state - weighted_u_mean)^2 *
        matrix(pi_state, nrow(draws$u_state), length(pi_state), byrow = TRUE)
    ))
    p_a <- draws$p_state_A
    weighted_p_mean <- as.vector(p_a %*% pi_state)
    sd_probability_w <- sqrt(rowSums(
      (p_a - weighted_p_mean)^2 *
        matrix(pi_state, nrow(p_a), length(pi_state), byrow = TRUE)
    ))
    .bhf_expect_close(draws$sd_state_logit_unweighted,
                      apply(draws$u_state, 1, stats::sd),
                      "unweighted logit finite-population SD")
    .bhf_expect_close(draws$sd_state_logit_weighted, sd_logit_w,
                      "weighted logit finite-population SD")
    .bhf_expect_close(draws$sd_state_probability_A_unweighted,
                      apply(p_a, 1, stats::sd),
                      "unweighted A probability finite-population SD")
    .bhf_expect_close(draws$sd_state_probability_A_weighted,
                      sd_probability_w,
                      "weighted A probability finite-population SD")

    raw_log_lik <- matrix(
      stats::dbinom(
        as.vector(matrix(prepared$stan_data$y,
                         nrow(draws$eta), prepared$stan_data$N, byrow = TRUE)),
        size = 1,
        prob = as.vector(stats::plogis(draws$eta)),
        log = TRUE
      ),
      nrow = nrow(draws$eta)
    )
    .bhf_expect_close(draws$log_lik_raw, raw_log_lik,
                      "raw pointwise log-likelihood")
    .bhf_expect_close(
      draws$log_lik_pseudo,
      raw_log_lik * matrix(prepared$stan_data$w_lik,
                           nrow(raw_log_lik), ncol(raw_log_lik), byrow = TRUE),
      "weighted pointwise pseudo-log-likelihood"
    )
    expect_true(all(draws$y_rep %in% c(0L, 1L)))
    expect_true(all(vapply(draws, function(x) all(is.finite(x)), logical(1))))

    make_parameter_draw <- function(alpha, z_state, z_stratum, z_psu,
                                    sigma_state, sigma_stratum, sigma_psu) {
      values <- c(
        alpha = alpha,
        stats::setNames(z_state, paste0("z_state[", seq_along(z_state), "]")),
        stats::setNames(
          z_stratum, paste0("z_stratum[", seq_along(z_stratum), "]")
        ),
        stats::setNames(z_psu, paste0("z_psu[", seq_along(z_psu), "]")),
        sigma_state = sigma_state,
        sigma_stratum = sigma_stratum,
        sigma_psu = sigma_psu
      )
      matrix(values, nrow = 1L, dimnames = list("edge", names(values)))
    }
    run_gqs <- function(data, parameter_draw, seed) {
      gqs_fit <- rstan::gqs(
        model, data = data, draws = parameter_draw, seed = seed
      )
      rstan::extract(gqs_fit, permuted = TRUE)
    }

    edge_data <- prepared$stan_data
    edge_data$w_state_pop_share <- c(0.5, 0.25, 0.25)
    edge_parameter <- make_parameter_draw(
      alpha = -0.4,
      z_state = c(1.2, -0.7, 0.3),
      z_stratum = c(-0.5, 0.5),
      z_psu = c(0.8, -0.8, 0.2, -0.2),
      sigma_state = 0.7,
      sigma_stratum = 0.3,
      sigma_psu = 0.4
    )

    disabled_data <- edge_data
    disabled_data$use_deattenuation <- 0L
    disabled_data$vhat_state <- rep(0, disabled_data$S)
    disabled <- run_gqs(disabled_data, edge_parameter, 520515)
    .bhf_expect_close(disabled$var_between_A_star, disabled$var_between_A,
                      "disabled A-star between carrier", tolerance = 0)
    .bhf_expect_close(disabled$var_within_A_star, disabled$var_within_A,
                      "disabled A-star within carrier", tolerance = 0)
    .bhf_expect_close(disabled$var_total_A_star, disabled$var_total_A,
                      "disabled A-star total carrier", tolerance = 0)
    .bhf_expect_close(disabled$prop_between_A_star, disabled$prop_between_A,
                      "disabled A-star proportion carrier", tolerance = 0)
    expect_identical(as.integer(disabled$A_star_at_boundary), 0L)
    expect_identical(as.integer(disabled$A_star_truncated), 0L)

    zero_vhat_data <- edge_data
    zero_vhat_data$use_deattenuation <- 1L
    zero_vhat_data$vhat_state <- rep(0, zero_vhat_data$S)
    zero_vhat <- run_gqs(zero_vhat_data, edge_parameter, 520516)
    boundary_data <- zero_vhat_data
    boundary_data$vhat_state <- rep(
      as.numeric(zero_vhat$var_between_A), boundary_data$S
    )
    boundary <- run_gqs(boundary_data, edge_parameter, 520517)
    expect_identical(as.numeric(boundary$var_between_A_star), 0)
    expect_identical(as.integer(boundary$A_star_at_boundary), 1L)
    expect_identical(as.integer(boundary$A_star_truncated), 0L)

    overshoot_data <- boundary_data
    overshoot_data$vhat_state <-
      boundary_data$vhat_state + rep(0.01, overshoot_data$S)
    overshoot <- run_gqs(overshoot_data, edge_parameter, 520518)
    expect_identical(as.numeric(overshoot$var_between_A_star), 0)
    expect_identical(as.integer(overshoot$A_star_at_boundary), 1L)
    expect_identical(as.integer(overshoot$A_star_truncated), 1L)

    zero_parameter <- make_parameter_draw(
      alpha = -1000,
      z_state = rep(0, edge_data$S),
      z_stratum = rep(0, edge_data$H),
      z_psu = rep(0, edge_data$J),
      sigma_state = 0,
      sigma_stratum = 0,
      sigma_psu = 0
    )
    total_zero_data <- disabled_data
    total_zero <- run_gqs(total_zero_data, zero_parameter, 520519)
    expect_identical(as.numeric(total_zero$var_total_A), 0)
    expect_identical(as.numeric(total_zero$prop_between_A), 0)
    expect_identical(as.integer(total_zero$A_proportion_defined), 0L)
    expect_identical(as.numeric(total_zero$var_total_B), 0)
    expect_identical(as.numeric(total_zero$prop_between_B), 0)
    expect_identical(as.integer(total_zero$B_proportion_defined), 0L)

    singleton_data <- list(
      data_schema_version = 1L,
      N = 2L, S = 1L, H = 1L, J = 1L,
      y = as.array(c(0L, 1L)),
      state_id = as.array(c(1L, 1L)),
      stratum_id = as.array(c(1L, 1L)),
      J_h = as.array(1L),
      psu_start = as.array(1L),
      psu_in_stratum_id = as.array(c(1L, 1L)),
      psu_flat_id = as.array(c(1L, 1L)),
      w_lik = as.array(c(1, 1)),
      w_state_pop_share = as.array(1),
      prior_alpha_mean = 0,
      prior_alpha_sd = 0.5,
      sigma_state_prior_code = 1L,
      use_deattenuation = 0L,
      vhat_state = as.array(0)
    )
    singleton_parameter <- make_parameter_draw(
      alpha = 0,
      z_state = 2,
      z_stratum = 1,
      z_psu = 3,
      sigma_state = 0.7,
      sigma_stratum = 0.4,
      sigma_psu = 0.2
    )
    singleton <- run_gqs(singleton_data, singleton_parameter, 520520)
    expect_identical(as.numeric(singleton$sd_state_logit_unweighted), 0)
    expect_identical(as.numeric(singleton$sd_state_logit_weighted), 0)
    expect_identical(
      as.numeric(singleton$sd_state_probability_A_unweighted), 0
    )
    expect_identical(
      as.numeric(singleton$sd_state_probability_A_weighted), 0
    )
    expect_identical(as.numeric(singleton$state_weight_B), 2)
  }
})
