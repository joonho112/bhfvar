// Bayesian Hybrid Framework for Variance Decomposition in Complex Surveys
// 
// This Stan model implements the methodology described in:
// Lee, J., & Hooper, A. (2025). Disentangling Signal from Noise: A Bayesian 
// Hybrid Framework for Variance Decomposition in Complex Surveys with Post-Hoc Domains.
//
// Key Features:
// - Bayesian Pseudo-Likelihood for design consistency
// - Hybrid GLMM with substantive (state) and nuisance (PSU) random effects
// - Dual Estimand Framework: Policy (A/A*) and Descriptive (B) estimands
// - De-attenuation for finite-sample variance inflation correction
//
// Author: JoonHo Lee (jlee296@ua.edu)
// Date: 2025-01-24

data {
  // Dimensions
  int<lower=1> N;                           // Total sample size
  int<lower=1> S;                           // Number of domains (states)
  int<lower=1> H;                           // Number of strata
  int<lower=1> J;                           // Total number of PSUs
  
  // Outcome (binary)
  array[N] int<lower=0, upper=1> y;         // Binary outcome
  
  // Hierarchical structure indices (1-indexed)
  array[N] int<lower=1, upper=S> state_id;  // State/domain indicator
  array[N] int<lower=1, upper=H> stratum_id; // Stratum indicator
  
  // PSU structure within strata
  array[H] int<lower=1> J_h;                // Number of PSUs in each stratum
  array[H] int<lower=1> psu_start;          // Starting index for PSUs in each stratum
  array[N] int<lower=1> psu_in_stratum_id;  // PSU index within stratum (1 to J_h[h])
  
  // Survey weights (scaled for pseudo-likelihood)
  array[N] real<lower=0> w_lik;             // Likelihood weights (scaled)
  
  // Population shares for states (sum to 1)
  array[S] real<lower=0, upper=1> w_state_pop_share;
  
  // Priors
  real prior_alpha_mean;                    // Prior mean for intercept
  real<lower=0> prior_alpha_sd;             // Prior SD for intercept
  
  // De-attenuation settings
  int<lower=0, upper=1> use_deattenuation;  // Whether to apply de-attenuation
  array[S] real<lower=0> vhat_state;        // Estimated sampling variance per state
}

transformed data {
  // Precompute normalized weights for pseudo-likelihood
  array[N] real w_norm;
  real w_sum = 0;
  
  for (i in 1:N) {
    w_sum += w_lik[i];
  }
  
  for (i in 1:N) {
    w_norm[i] = w_lik[i] / w_sum * N;
  }
  
  // Compute population-weighted mean of vhat for de-attenuation
  real vhat_mean = 0;
  if (use_deattenuation == 1) {
    for (s in 1:S) {
      vhat_mean += w_state_pop_share[s] * vhat_state[s];
    }
  }
}

parameters {
  // Fixed effects
  real alpha;                               // Global intercept (logit scale)
  
  // Random effects (non-centered parameterization)
  vector[S] z_state;                        // Standardized state effects
  vector[J] z_psu;                          // Standardized PSU effects
  
  // Variance components
  real<lower=0> sigma_state;                // Between-state SD (logit scale)
  real<lower=0> sigma_psu;                  // Between-PSU SD (within stratum)
}

transformed parameters {
  // State random effects
  vector[S] u_state = sigma_state * z_state;
  
  // PSU random effects  
  vector[J] u_psu = sigma_psu * z_psu;
  
  // Linear predictor
  vector[N] eta;
  for (i in 1:N) {
    int h = stratum_id[i];
    int psu_idx = psu_start[h] + psu_in_stratum_id[i] - 1;
    eta[i] = alpha + u_state[state_id[i]] + u_psu[psu_idx];
  }
  
  // ============================================================
  // Variance Components on Logit Scale (Estimand A - Policy)
  // ============================================================
  real var_between_state = square(sigma_state);
  real var_psu = square(sigma_psu);
  real var_logistic = square(pi()) / 3;  // Level-1 variance for logistic
  real var_within_state = var_psu + var_logistic;
  real var_total_logit = var_between_state + var_within_state;
  
  // ICC on logit scale (Estimand A)
  real icc_state = var_between_state / var_total_logit;
  
  // Marginal scaling factor (Zeger et al., 1988)
  real c_factor = sqrt(1.0 / (1.0 + 0.346 * var_between_state));
}

model {
  // ============================================================
  // Priors
  // ============================================================
  
  // Intercept prior (informative based on overall proportion)
  alpha ~ normal(prior_alpha_mean, prior_alpha_sd);
  
  // Variance component priors (half-normal)
  sigma_state ~ normal(0, 1);
  sigma_psu ~ normal(0, 0.5);
  
  // Standard normal priors for non-centered effects
  z_state ~ std_normal();
  z_psu ~ std_normal();
  
  // ============================================================
  // Pseudo-Likelihood
  // ============================================================
  for (i in 1:N) {
    target += w_norm[i] * bernoulli_logit_lpmf(y[i] | eta[i]);
  }
}

generated quantities {
  // ============================================================
  // State-Level Probabilities
  // ============================================================
  
  // Conditional probabilities (given state effect)
  array[S] real p_state_conditional;
  for (s in 1:S) {
    p_state_conditional[s] = inv_logit(alpha + u_state[s]);
  }
  
  // Marginal probabilities (integrating out within-state variation)
  array[S] real p_state_marginal;
  for (s in 1:S) {
    real eta_marginal = (alpha + u_state[s]) * c_factor;
    p_state_marginal[s] = inv_logit(eta_marginal);
  }
  
  // Population-weighted overall probability
  real p_overall = 0;
  for (s in 1:S) {
    p_overall += w_state_pop_share[s] * p_state_marginal[s];
  }
  
  // ============================================================
  // Variance Decomposition on Probability Scale (Estimand B - Descriptive)
  // ============================================================
  
  // Between-state variance (finite population)
  real var_between_prob = 0;
  for (s in 1:S) {
    var_between_prob += w_state_pop_share[s] * square(p_state_marginal[s] - p_overall);
  }
  
  // Within-state variance (Bernoulli variance)
  real var_within_prob = 0;
  for (s in 1:S) {
    var_within_prob += w_state_pop_share[s] * p_state_marginal[s] * (1 - p_state_marginal[s]);
  }
  
  // Total variance
  real var_total_prob = var_between_prob + var_within_prob;
  
  // ICC on probability scale (Estimand B)
  real icc_prob = var_between_prob / var_total_prob;
  
  // ============================================================
  // De-attenuated Estimands (Estimand A* - Policy Adjusted)
  // ============================================================
  
  // De-attenuated between-state variance
  real var_between_deatten;
  real icc_deatten;
  
  if (use_deattenuation == 1) {
    // Subtract estimated sampling variance (with floor at small positive value)
    var_between_deatten = fmax(0.001, var_between_prob - vhat_mean);
    
    // Recompute ICC with de-attenuated variance
    real var_total_deatten = var_between_deatten + var_within_prob;
    icc_deatten = var_between_deatten / var_total_deatten;
  } else {
    var_between_deatten = var_between_prob;
    icc_deatten = icc_prob;
  }
  
  // ============================================================
  // Reliability and Effective Sample Size
  // ============================================================
  
  // State-specific reliability (shrinkage factor)
  array[S] real reliability_state;
  for (s in 1:S) {
    // Using the formula: R = sigma^2 / (sigma^2 + V_s)
    // where V_s is the sampling variance for state s
    if (vhat_state[s] > 0) {
      reliability_state[s] = var_between_prob / (var_between_prob + vhat_state[s]);
    } else {
      reliability_state[s] = 1.0;
    }
  }
  
  // Average reliability (population-weighted)
  real reliability_avg = 0;
  for (s in 1:S) {
    reliability_avg += w_state_pop_share[s] * reliability_state[s];
  }
  
  // ============================================================
  // Model Fit Diagnostics
  // ============================================================
  
  // Log-likelihood for LOO-CV
  array[N] real log_lik;
  for (i in 1:N) {
    log_lik[i] = bernoulli_logit_lpmf(y[i] | eta[i]);
  }
  
  // Posterior predictive
  array[N] int y_rep;
  for (i in 1:N) {
    y_rep[i] = bernoulli_logit_rng(eta[i]);
  }
}
