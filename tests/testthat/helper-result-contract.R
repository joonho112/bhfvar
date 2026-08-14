make_result_contract_draws <- function(enabled = TRUE) {
  D <- 5L
  share <- c(0.5, 0.3, 0.2)
  p_state_A <- rbind(
    c(0.20, 0.40, 0.60), c(0.25, 0.45, 0.65),
    c(0.30, 0.50, 0.70), c(0.35, 0.55, 0.75),
    c(0.40, 0.60, 0.80)
  )
  p_state_B <- p_state_A + matrix(
    c(0.02, -0.01, 0.03), D, 3L, byrow = TRUE
  )
  p_bar_A <- as.vector(p_state_A %*% share)
  p_bar_B <- as.vector(p_state_B %*% share)
  between_A <- vapply(seq_len(D), function(draw) {
    sum(share * (p_state_A[draw, ] - p_bar_A[draw])^2)
  }, numeric(1))
  within_A <- vapply(seq_len(D), function(draw) {
    sum(share * p_state_A[draw, ] * (1 - p_state_A[draw, ]))
  }, numeric(1))
  total_A <- between_A + within_A
  proportion_A <- between_A / total_A

  within_binomial_B <- seq(0.16, 0.14, length.out = D)
  within_mixture_B <- seq(0.01, 0.014, length.out = D)
  within_B <- within_binomial_B + within_mixture_B
  between_B <- vapply(seq_len(D), function(draw) {
    sum(share * (p_state_B[draw, ] - p_bar_B[draw])^2)
  }, numeric(1))
  total_B <- between_B + within_B
  proportion_B <- between_B / total_B

  correction <- if (enabled) rep(0.005, D) else rep(0, D)
  between_A_star <- if (enabled) pmax(0, between_A - correction) else between_A
  total_A_star <- between_A_star + within_A
  proportion_A_star <- between_A_star / total_A_star

  raw_log_lik <- -matrix(seq_len(D * 6L) / 10, D, 6L)
  w_lik <- c(0.5, 1.5, 0.75, 1.25, 0.8, 1.2)
  pseudo_log_lik <- sweep(raw_log_lik, 2, w_lik, `*`)

  list(
    sigma_state = seq(0.4, 0.8, length.out = D),
    sigma_stratum = seq(0.3, 0.5, length.out = D),
    sigma_psu = seq(0.6, 1.0, length.out = D),
    var_state_latent = seq(0.16, 0.64, length.out = D),
    var_stratum_latent = seq(0.09, 0.25, length.out = D),
    var_psu_latent = seq(0.36, 1.00, length.out = D),
    var_level1_latent = rep(pi^2 / 3, D),
    var_total_latent = seq(0.16, 0.64, length.out = D) +
      seq(0.09, 0.25, length.out = D) +
      seq(0.36, 1.00, length.out = D) + pi^2 / 3,
    icc_state_latent = seq(0.05, 0.12, length.out = D),
    icc_stratum_latent = seq(0.03, 0.06, length.out = D),
    icc_psu_latent = seq(0.10, 0.18, length.out = D),
    p_state_A = p_state_A,
    p_bar_A = p_bar_A,
    var_between_A = between_A,
    var_within_A = within_A,
    var_total_A = total_A,
    prop_between_A = proportion_A,
    A_proportion_defined = rep(1L, D),
    p_bar_A_star = p_bar_A,
    var_between_A_star = between_A_star,
    var_within_A_star = within_A,
    var_total_A_star = total_A_star,
    prop_between_A_star = proportion_A_star,
    mean_vhat_state = correction,
    A_star_at_boundary = as.integer(between_A_star == 0),
    A_star_truncated = as.integer(enabled & between_A < correction),
    A_star_proportion_defined = rep(1L, D),
    p_state_B = p_state_B,
    p_bar_B = p_bar_B,
    var_between_B = between_B,
    var_within_binomial_B = within_binomial_B,
    var_within_mixture_B = within_mixture_B,
    var_within_B = within_B,
    var_total_B = total_B,
    prop_between_B = proportion_B,
    B_proportion_defined = rep(1L, D),
    gap_B_minus_A_mean = p_bar_B - p_bar_A,
    gap_B_minus_A_between = between_B - between_A,
    gap_B_minus_A_within = within_B - within_A,
    gap_B_minus_A_total = total_B - total_A,
    gap_B_minus_A_proportion = proportion_B - proportion_A,
    gap_A_minus_A_star_between = between_A - between_A_star,
    gap_A_minus_A_star_total = total_A - total_A_star,
    gap_A_minus_A_star_proportion = proportion_A - proportion_A_star,
    log_lik_raw = raw_log_lik,
    log_lik_pseudo = pseudo_log_lik
  )
}

make_result_contract_fit <- function(enabled = TRUE) {
  analysis_data <- data.frame(
    original_row = c(2L, 4L, 5L, 8L, 9L, 10L),
    state_id = c(1L, 1L, 2L, 2L, 3L, 3L),
    stratum_id = c(1L, 1L, 1L, 2L, 2L, 2L),
    psu_flat_id = c(1L, 1L, 2L, 3L, 4L, 4L)
  )
  data <- list(
    schema_version = "0.5.0",
    contract_id = "bhfvar-data-contract-0.5.0",
    stan_data_schema_version = 1L,
    stan_data = list(
      data_schema_version = 1L,
      N = 6L, S = 3L, H = 2L, J = 4L,
      state_id = as.array(analysis_data$state_id),
      stratum_id = as.array(analysis_data$stratum_id),
      psu_flat_id = as.array(analysis_data$psu_flat_id),
      w_lik = as.array(c(0.5, 1.5, 0.75, 1.25, 0.8, 1.2)),
      w_state_pop_share = as.array(c(0.5, 0.3, 0.2)),
      use_deattenuation = as.integer(enabled)
    ),
    mapping = list(
      domain = data.frame(label = c("A", "B", "C"), id = 1:3),
      stratum = data.frame(label = c("H1", "H2"), id = 1:2),
      psu = data.frame(
        label = c("P1", "P2", "P1", "P2"), id = 1:4,
        stratum_label = c("H1", "H1", "H2", "H2")
      )
    ),
    analysis_data = analysis_data,
    domain_summary = data.frame(state_id = 1:3, n = c(2L, 2L, 2L)),
    provenance = list(sampling_variances = list(
      mode = if (enabled) "supplied" else "none",
      fixed_input = enabled
    ))
  )
  class(data) <- c("bhf_data", "list")
  structure(
    list(
      schema_version = "0.5.0",
      contract_id = "bhfvar-fit-contract-0.5.0",
      stanfit = structure(list(), class = "mock_stanfit"),
      data = data,
      .draws = make_result_contract_draws(enabled)
    ),
    class = c("bhf_fit", "list")
  )
}

local_result_draw_backend <- function() {
  testthat::local_mocked_bindings(
    .bhf_extract_draws = function(fit, pars) fit$.draws[pars],
    .package = "bhfvar",
    .env = parent.frame()
  )
}
