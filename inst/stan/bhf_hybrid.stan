// Bayesian Hybrid Framework for Variance Decomposition in Complex Surveys
// Article-aligned core for Mathematics 14(3), 512 (2026).
// State is crossed with stratum; PSU is strictly nested within stratum.

data {
  int<lower=1> data_schema_version;
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> H;
  int<lower=1> J;

  array[N] int<lower=0, upper=1> y;
  array[N] int<lower=1, upper=S> state_id;
  array[N] int<lower=1, upper=H> stratum_id;

  array[H] int<lower=1> J_h;
  array[H] int<lower=1> psu_start;
  array[N] int<lower=1> psu_in_stratum_id;
  array[N] int<lower=1, upper=J> psu_flat_id;

  vector<lower=0>[N] w_lik;
  vector<lower=0, upper=1>[S] w_state_pop_share;

  real prior_alpha_mean;
  real<lower=0> prior_alpha_sd;
  int<lower=1, upper=4> sigma_state_prior_code;

  int<lower=0, upper=1> use_deattenuation;
  vector<lower=0>[S] vhat_state;
}

transformed data {
  {
    vector[S] state_weight_check = rep_vector(0.0, S);
    for (i in 1:N) {
      state_weight_check[state_id[i]] += w_lik[i];
    }
    if (min(state_weight_check) <= 0) {
      reject("every state must have positive likelihood-weight support");
    }
  }
  if (data_schema_version != 1) {
    reject("Unsupported data_schema_version: ", data_schema_version);
  }
  if (min(w_lik) <= 0) {
    reject("w_lik must be strictly positive");
  }
  if (abs(sum(w_lik) - N) > 1e-10 * fmax(1.0, N)) {
    reject("w_lik must use the canonical global sum-N normalization");
  }
  if (min(w_state_pop_share) <= 0) {
    reject("w_state_pop_share must be strictly positive");
  }
  if (abs(sum(w_state_pop_share) - 1.0) > 1e-12) {
    reject("w_state_pop_share must sum to one");
  }
  if (prior_alpha_sd <= 0) {
    reject("prior_alpha_sd must be strictly positive");
  }
  if (use_deattenuation == 0 && max(vhat_state) > 0) {
    reject("disabled de-attenuation requires the approved zero placeholder");
  }
  if (sum(J_h) != J) {
    reject("J_h must sum to J");
  }
  if (psu_start[1] != 1) {
    reject("the first PSU block must start at one");
  }
  if (H > 1) {
    for (h in 2:H) {
      if (psu_start[h] != psu_start[h - 1] + J_h[h - 1]) {
        reject("psu_start must define consecutive stratum blocks");
      }
    }
  }
  for (i in 1:N) {
    int h = stratum_id[i];
    if (psu_in_stratum_id[i] > J_h[h]) {
      reject("psu_in_stratum_id is outside its stratum block");
    }
    if (psu_flat_id[i] != psu_start[h] + psu_in_stratum_id[i] - 1) {
      reject("psu_flat_id does not match its stratum block reconstruction");
    }
  }
}

parameters {
  real alpha;

  vector[S] z_state;
  vector[H] z_stratum;
  vector[J] z_psu;

  real<lower=0> sigma_state;
  real<lower=0> sigma_stratum;
  real<lower=0> sigma_psu;
}

transformed parameters {
  vector[S] u_state;
  vector[H] u_stratum;
  vector[J] u_psu;
  vector[N] eta;

  {
    real state_weighted_mean = dot_product(w_state_pop_share, z_state);
    u_state = sigma_state *
      (z_state - rep_vector(state_weighted_mean, S));
  }

  u_stratum = sigma_stratum *
    (z_stratum - rep_vector(mean(z_stratum), H));

  for (h in 1:H) {
    real block_mean = 0;
    for (j_local in 1:J_h[h]) {
      int j_flat = psu_start[h] + j_local - 1;
      block_mean += z_psu[j_flat];
    }
    block_mean /= J_h[h];
    for (j_local in 1:J_h[h]) {
      int j_flat = psu_start[h] + j_local - 1;
      u_psu[j_flat] = sigma_psu * (z_psu[j_flat] - block_mean);
    }
  }

  for (i in 1:N) {
    eta[i] = alpha +
      u_state[state_id[i]] +
      u_stratum[stratum_id[i]] +
      u_psu[psu_flat_id[i]];
  }
}

model {
  alpha ~ normal(prior_alpha_mean, prior_alpha_sd);

  z_state ~ std_normal();
  z_stratum ~ std_normal();
  z_psu ~ std_normal();

  if (sigma_state_prior_code == 1) {
    sigma_state ~ student_t(3, 0, 2.5);
  } else if (sigma_state_prior_code == 2) {
    sigma_state ~ normal(0, 1);
  } else if (sigma_state_prior_code == 3) {
    sigma_state ~ cauchy(0, 2.5);
  } else {
    sigma_state ~ student_t(3, 0, 5);
  }
  sigma_stratum ~ student_t(3, 0, 2.5);
  sigma_psu ~ student_t(3, 0, 2.5);

  for (i in 1:N) {
    target += w_lik[i] * bernoulli_logit_lpmf(y[i] | eta[i]);
  }
}

generated quantities {
  // Phase 4 centering diagnostics retained as a structural debug surface.
  real state_centering_residual;
  real stratum_centering_residual;
  vector[H] psu_centering_residual;

  // Latent-scale diagnostics. These are not probability-scale estimands.
  real var_state_latent;
  real var_stratum_latent;
  real var_psu_latent;
  real var_level1_latent;
  real var_total_latent;
  real icc_state_latent;
  real icc_stratum_latent;
  real icc_psu_latent;

  // Estimand A: state-only policy decomposition.
  vector[S] p_state_A;
  real p_bar_A;
  real var_between_A;
  real var_within_A;
  real var_total_A;
  real prop_between_A;
  int<lower=0, upper=1> A_proportion_defined;

  // Estimand A*: de-attenuated A. When disabled, these finite carriers equal A;
  // the R schema marks A* unavailable and does not present them as computed.
  real mean_vhat_state;
  real p_bar_A_star;
  real var_between_A_star;
  real var_within_A_star;
  real var_total_A_star;
  real prop_between_A_star;
  int<lower=0, upper=1> A_star_at_boundary;
  int<lower=0, upper=1> A_star_truncated;
  int<lower=0, upper=1> A_star_proportion_defined;

  // Estimand B: realized full-design descriptive decomposition.
  vector[N] p_individual_B;
  vector[S] state_weight_B;
  vector[S] p_state_B;
  real p_bar_B;
  vector[S] within_binomial_state_B;
  vector[S] within_mixture_state_B;
  real var_within_binomial_B;
  real var_within_mixture_B;
  real var_within_B;
  real var_between_B;
  real var_total_B;
  real prop_between_B;
  int<lower=0, upper=1> B_proportion_defined;

  // Signed diagnostic gaps.
  real gap_B_minus_A_mean;
  real gap_B_minus_A_between;
  real gap_B_minus_A_within;
  real gap_B_minus_A_total;
  real gap_B_minus_A_proportion;
  real gap_A_minus_A_star_between;
  real gap_A_minus_A_star_total;
  real gap_A_minus_A_star_proportion;

  // Appendix finite-population state dispersion diagnostics.
  real sd_state_logit_unweighted;
  real sd_state_logit_weighted;
  real sd_state_probability_A_unweighted;
  real sd_state_probability_A_weighted;

  // Pointwise raw and pseudo log likelihoods plus unweighted replication.
  vector[N] log_lik_pseudo;
  vector[N] log_lik_raw;
  array[N] int<lower=0, upper=1> y_rep;

  state_centering_residual = dot_product(w_state_pop_share, u_state);
  stratum_centering_residual = mean(u_stratum);
  for (h in 1:H) {
    real block_sum = 0;
    for (j_local in 1:J_h[h]) {
      int j_flat = psu_start[h] + j_local - 1;
      block_sum += u_psu[j_flat];
    }
    psu_centering_residual[h] = block_sum / J_h[h];
  }

  var_state_latent = square(sigma_state);
  var_stratum_latent = square(sigma_stratum);
  var_psu_latent = square(sigma_psu);
  var_level1_latent = square(pi()) / 3.0;
  var_total_latent = var_state_latent + var_stratum_latent +
    var_psu_latent + var_level1_latent;
  icc_state_latent = var_state_latent / var_total_latent;
  icc_stratum_latent = var_stratum_latent / var_total_latent;
  icc_psu_latent = var_psu_latent / var_total_latent;

  for (s in 1:S) {
    p_state_A[s] = inv_logit(alpha + u_state[s]);
  }
  p_bar_A = dot_product(w_state_pop_share, p_state_A);
  var_between_A = 0.0;
  var_within_A = 0.0;
  for (s in 1:S) {
    var_between_A += w_state_pop_share[s] *
      square(p_state_A[s] - p_bar_A);
    var_within_A += w_state_pop_share[s] * p_state_A[s] *
      (1.0 - p_state_A[s]);
  }
  var_total_A = var_between_A + var_within_A;
  A_proportion_defined = var_total_A > 0;
  prop_between_A = A_proportion_defined == 1
    ? var_between_A / var_total_A
    : 0.0;

  if (use_deattenuation == 1) {
    real raw_between_A_star;
    mean_vhat_state = dot_product(w_state_pop_share, vhat_state);
    raw_between_A_star = var_between_A - mean_vhat_state;
    var_between_A_star = fmax(0.0, raw_between_A_star);
    A_star_at_boundary = var_between_A_star == 0.0;
    A_star_truncated = raw_between_A_star < 0.0;
  } else {
    mean_vhat_state = 0.0;
    var_between_A_star = var_between_A;
    A_star_at_boundary = 0;
    A_star_truncated = 0;
  }
  p_bar_A_star = p_bar_A;
  var_within_A_star = var_within_A;
  var_total_A_star = var_between_A_star + var_within_A_star;
  A_star_proportion_defined = var_total_A_star > 0;
  prop_between_A_star = A_star_proportion_defined == 1
    ? var_between_A_star / var_total_A_star
    : 0.0;

  // First observation pass: probabilities, state weights, state probability
  // numerators, and Bernoulli numerators.
  state_weight_B = rep_vector(0.0, S);
  p_state_B = rep_vector(0.0, S);
  within_binomial_state_B = rep_vector(0.0, S);
  within_mixture_state_B = rep_vector(0.0, S);
  for (i in 1:N) {
    int s = state_id[i];
    real p = inv_logit(eta[i]);
    p_individual_B[i] = p;
    state_weight_B[s] += w_lik[i];
    p_state_B[s] += w_lik[i] * p;
    within_binomial_state_B[s] += w_lik[i] * p * (1.0 - p);
  }
  for (s in 1:S) {
    p_state_B[s] /= state_weight_B[s];
    within_binomial_state_B[s] /= state_weight_B[s];
  }

  // Second observation pass: deviation-form mixture variance. This avoids the
  // cancellation-prone E[p^2] - E[p]^2 form while remaining O(N + S).
  for (i in 1:N) {
    int s = state_id[i];
    within_mixture_state_B[s] += w_lik[i] *
      square(p_individual_B[i] - p_state_B[s]);
  }
  for (s in 1:S) {
    within_mixture_state_B[s] /= state_weight_B[s];
  }

  p_bar_B = dot_product(w_state_pop_share, p_state_B);
  var_within_binomial_B = dot_product(
    w_state_pop_share, within_binomial_state_B
  );
  var_within_mixture_B = dot_product(
    w_state_pop_share, within_mixture_state_B
  );
  var_within_B = var_within_binomial_B + var_within_mixture_B;
  var_between_B = 0.0;
  for (s in 1:S) {
    var_between_B += w_state_pop_share[s] *
      square(p_state_B[s] - p_bar_B);
  }
  var_total_B = var_between_B + var_within_B;
  B_proportion_defined = var_total_B > 0;
  prop_between_B = B_proportion_defined == 1
    ? var_between_B / var_total_B
    : 0.0;

  gap_B_minus_A_mean = p_bar_B - p_bar_A;
  gap_B_minus_A_between = var_between_B - var_between_A;
  gap_B_minus_A_within = var_within_B - var_within_A;
  gap_B_minus_A_total = var_total_B - var_total_A;
  gap_B_minus_A_proportion = prop_between_B - prop_between_A;
  gap_A_minus_A_star_between = var_between_A - var_between_A_star;
  gap_A_minus_A_star_total = var_total_A - var_total_A_star;
  gap_A_minus_A_star_proportion =
    prop_between_A - prop_between_A_star;

  sd_state_logit_unweighted = S > 1 ? sd(u_state) : 0.0;
  sd_state_probability_A_unweighted = S > 1 ? sd(p_state_A) : 0.0;
  {
    real mean_logit_weighted = dot_product(w_state_pop_share, u_state);
    real var_logit_weighted = 0.0;
    real var_probability_A_weighted = 0.0;
    for (s in 1:S) {
      var_logit_weighted += w_state_pop_share[s] *
        square(u_state[s] - mean_logit_weighted);
      var_probability_A_weighted += w_state_pop_share[s] *
        square(p_state_A[s] - p_bar_A);
    }
    sd_state_logit_weighted = sqrt(fmax(0.0, var_logit_weighted));
    sd_state_probability_A_weighted =
      sqrt(fmax(0.0, var_probability_A_weighted));
  }

  for (i in 1:N) {
    log_lik_raw[i] = bernoulli_logit_lpmf(y[i] | eta[i]);
    log_lik_pseudo[i] = w_lik[i] * log_lik_raw[i];
    y_rep[i] = bernoulli_logit_rng(eta[i]);
  }
}
