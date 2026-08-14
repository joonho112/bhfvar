# Phase 7 scientific-validation plans and collectors.

.bhf_validation_hash <- function(x) {
  if (exists(".bhf_synthetic_checksum", mode = "function")) {
    return(.bhf_synthetic_checksum(x))
  }
  text <- paste(utils::capture.output(base::dput(x)), collapse = "\n")
  path <- tempfile("bhf-validation-hash-")
  on.exit(unlink(path), add = TRUE)
  writeLines(text, path, useBytes = TRUE)
  unname(as.character(tools::md5sum(path)))
}

#' Frozen Scientific Validation Plan
#'
#' Internal plan mirroring Gate G1. It is data/config only and launches no job.
#' @keywords internal
bhf_scientific_validation_plan <- function() {
  medium <- expand.grid(
    profile = c("low", "high"),
    rho = c(0, 0.5),
    seed_index = seq_len(3L),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  medium$effect_seed <- 720000L +
    match(medium$profile, c("low", "high")) * 1000L +
    as.integer(medium$rho * 100L) + medium$seed_index * 10L
  medium$outcome_seed <- medium$effect_seed + 1L
  medium$mcmc_seed <- medium$effect_seed + 2L
  medium$job_id <- sprintf(
    "recovery-%s-rho%s-seed%d",
    medium$profile, ifelse(medium$rho == 0, "0", "05"), medium$seed_index
  )
  medium <- medium[c(
    "job_id", "profile", "rho", "seed_index", "effect_seed",
    "outcome_seed", "mcmc_seed"
  )]

  plan <- list(
    schema_version = "bhfvar-scientific-validation-plan-1.0.0",
    authority = "Gate G1 frozen stochastic criteria",
    tiny = list(
      job_id = "tiny-low-rho05-seed7109001",
      profile = "low", rho = 0.5,
      effect_seed = 7106001L, outcome_seed = 7107001L,
      mcmc_seed = 7109001L,
      dgp = list(n_states = 8L, n_strata = 4L,
                 psus_per_stratum = 2L, observations_per_cell = 2L),
      mcmc = list(chains = 4L, iter = 1000L, warmup = 500L,
                  kept_per_chain = 500L, adapt_delta = 0.95,
                  max_treedepth = 12L),
      thresholds = list(
        divergences = 0L, treedepth_hits = 0L, ebfmi_min = 0.20,
        rhat_max = 1.05, ess_bulk_min = 100, ess_tail_min = 100
      )
    ),
    medium = list(
      jobs = medium,
      dgp = list(n_states = 10L, n_strata = 6L,
                 psus_per_stratum = 2L, observations_per_cell = 2L),
      mcmc = list(chains = 4L, iter = 2000L, warmup = 1000L,
                  kept_per_chain = 1000L, adapt_delta = 0.95,
                  max_treedepth = 12L),
      thresholds = list(
        divergences = 0L, treedepth_hits = 0L, ebfmi_min = 0.20,
        core_rhat_max = 1.01, ess_bulk_min = 400,
        ess_tail_min = 400, raw_effect_rhat_max = 1.05,
        interval_probability = 0.90, overall_coverage_min = 0.80,
        target_family_coverage_min = 8L, standardized_error_median_max = 1.0,
        standardized_error_q95_max = 2.5,
        informative_gap_sign_probability_min = 0.80,
        informative_gap_seed_pass_min = 2L
      )
    ),
    constraints = list(
      oracle_immutable = TRUE,
      sensitivity_cannot_rescue_recovery = TRUE,
      restricted_application_reproduction_claim = FALSE,
      doubled_iteration_rerun_max = 1L
    )
  )
  plan$config_hash <- .bhf_synthetic_checksum(plan)
  plan
}

#' Scientific Sensitivity Grid
#' @keywords internal
bhf_scientific_sensitivity_grid <- function() {
  data.frame(
    variant = c(
      "half_t3_2.5", "half_normal_1", "half_cauchy_2.5", "half_t3_5"
    ),
    sigma_state_prior_code = 1:4,
    varied_component = rep("sigma_state", 4L),
    other_priors_fixed = rep(TRUE, 4L),
    stringsAsFactors = FALSE
  )
}

.bhf_validation_generator <- function() {
  if (exists("generate_bhf_synthetic_dgp", mode = "function")) {
    return(get("generate_bhf_synthetic_dgp", mode = "function"))
  }
  if (exists("generate_article_synthetic", mode = "function")) {
    return(get("generate_article_synthetic", mode = "function"))
  }
  stop("No article-aligned synthetic generator is available.", call. = FALSE)
}

#' Materialize a Validation Job Without Sampling
#' @keywords internal
bhf_materialize_validation_job <- function(job, dgp, prior = "half_t3_2.5") {
  required <- c("job_id", "profile", "rho", "effect_seed", "outcome_seed")
  if (!is.list(job) || !all(required %in% names(job))) {
    stop("job is missing a frozen validation field.", call. = FALSE)
  }
  generator <- .bhf_validation_generator()
  synthetic <- do.call(generator, c(
    list(
      profile = job$profile, rho = job$rho,
      effect_seed = as.integer(job$effect_seed),
      outcome_seed = as.integer(job$outcome_seed)
    ),
    dgp
  ))
  prepared <- prepare_bhf_data(
    synthetic$data,
    outcome = "outcome", domain = "state", strata = "stratum",
    psu = "psu", weights = "weight",
    population_shares = synthetic$truth$population_shares,
    deattenuation = "supplied",
    sampling_variances = synthetic$truth$vhat_state,
    sampling_variance_method = "external_other",
    prior_alpha_mean = synthetic$truth$config$alpha,
    prior_alpha_sd = 0.5,
    sigma_state_prior = prior
  )
  list(job = job, synthetic = synthetic, prepared = prepared)
}

.bhf_validation_ebfmi <- function(stanfit) {
  sampler <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)
  vapply(sampler, function(chain) {
    energy <- as.numeric(chain[, "energy__"])
    if (length(energy) < 2L || !all(is.finite(energy)) ||
        stats::var(energy) <= 0) return(NA_real_)
    mean(diff(energy)^2) / stats::var(energy)
  }, numeric(1))
}

.bhf_validation_draw_diagnostics <- function(stanfit, variables) {
  raw <- rstan::extract(
    stanfit, pars = variables, permuted = FALSE, inc_warmup = FALSE
  )
  draws <- posterior::as_draws_array(raw)
  summary <- posterior::summarise_draws(
    draws, "rhat", "ess_bulk", "ess_tail"
  )
  as.data.frame(summary, stringsAsFactors = FALSE)
}

#' Collect Frozen MCMC Diagnostics
#' @keywords internal
bhf_collect_validation_diagnostics <- function(fit, tier = c("tiny", "medium")) {
  tier <- match.arg(tier)
  if (!inherits(fit, "bhf_fit")) stop("fit must be a bhf_fit.", call. = FALSE)
  sampler <- rstan::get_sampler_params(fit$stanfit, inc_warmup = FALSE)
  divergences <- sum(vapply(
    sampler, function(x) sum(x[, "divergent__"]), numeric(1)
  ))
  depth <- fit$provenance$sampling$control$max_treedepth
  treedepth_hits <- sum(vapply(
    sampler, function(x) sum(x[, "treedepth__"] >= depth), numeric(1)
  ))
  core <- c("alpha", "sigma_state", "sigma_stratum", "sigma_psu")
  raw <- c("z_state", "z_stratum", "z_psu")
  list(
    tier = tier,
    divergences = as.integer(divergences),
    treedepth_hits = as.integer(treedepth_hits),
    ebfmi = .bhf_validation_ebfmi(fit$stanfit),
    core = .bhf_validation_draw_diagnostics(fit$stanfit, core),
    raw_effect = .bhf_validation_draw_diagnostics(fit$stanfit, raw)
  )
}

.bhf_validation_scalar_targets <- c(
  "alpha", "sigma_state", "sigma_stratum", "sigma_psu",
  "p_bar_A", "var_between_A", "var_within_A", "prop_between_A",
  "p_bar_A_star", "var_between_A_star", "var_within_A_star",
  "prop_between_A_star", "p_bar_B", "var_between_B", "var_within_B",
  "prop_between_B", "gap_B_minus_A_between"
)

#' Collect Scientific Target Draws
#' @keywords internal
bhf_collect_validation_targets <- function(fit) {
  if (!inherits(fit, "bhf_fit")) stop("fit must be a bhf_fit.", call. = FALSE)
  draws <- rstan::extract(
    fit$stanfit, pars = .bhf_validation_scalar_targets, permuted = TRUE
  )
  missing <- setdiff(.bhf_validation_scalar_targets, names(draws))
  if (length(missing)) {
    stop("fit is missing validation targets: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  as.data.frame(lapply(draws, as.numeric), check.names = FALSE)
}

#' Adjudicate Per-Fit Diagnostics Against Gate G1
#' @keywords internal
bhf_adjudicate_validation_diagnostics <- function(diagnostics,
                                                   tier = diagnostics$tier) {
  plan <- bhf_scientific_validation_plan()
  threshold <- if (identical(tier, "tiny")) {
    plan$tiny$thresholds
  } else {
    plan$medium$thresholds
  }
  core_rhat_limit <- if (identical(tier, "tiny")) threshold$rhat_max else
    threshold$core_rhat_max
  checks <- c(
    divergences = diagnostics$divergences == threshold$divergences,
    treedepth = diagnostics$treedepth_hits == threshold$treedepth_hits,
    ebfmi = all(is.finite(diagnostics$ebfmi)) &&
      min(diagnostics$ebfmi) >= threshold$ebfmi_min,
    core_rhat = length(diagnostics$core$rhat) > 0L &&
      all(is.finite(diagnostics$core$rhat)) &&
      max(diagnostics$core$rhat) <= core_rhat_limit,
    bulk_ess = length(diagnostics$core$ess_bulk) > 0L &&
      all(is.finite(diagnostics$core$ess_bulk)) &&
      min(diagnostics$core$ess_bulk) >= threshold$ess_bulk_min,
    tail_ess = length(diagnostics$core$ess_tail) > 0L &&
      all(is.finite(diagnostics$core$ess_tail)) &&
      min(diagnostics$core$ess_tail) >= threshold$ess_tail_min
  )
  if (identical(tier, "medium")) {
    checks <- c(
      checks,
      raw_effect_rhat = length(diagnostics$raw_effect$rhat) > 0L &&
        all(is.finite(diagnostics$raw_effect$rhat)) &&
        max(diagnostics$raw_effect$rhat) <= threshold$raw_effect_rhat_max
    )
  }
  list(pass = all(checks), checks = checks, thresholds = threshold)
}

#' Check Whether a Validation Artifact Can Be Resumed
#' @keywords internal
bhf_validation_artifact_status <- function(path, config_hash,
                                           expected = list()) {
  if (!file.exists(path)) return(list(status = "missing", reusable = FALSE))
  artifact <- tryCatch(readRDS(path), error = identity)
  if (inherits(artifact, "error") || !is.list(artifact)) {
    return(list(status = "corrupt", reusable = FALSE))
  }
  if (!identical(artifact$config_hash, config_hash) ||
      !identical(artifact$status, "complete")) {
    return(list(status = "stale_or_incomplete", reusable = FALSE))
  }
  if (!is.list(expected) ||
      (length(expected) &&
       (is.null(names(expected)) || any(!nzchar(names(expected)))))) {
    stop("expected must be an empty or fully named list", call. = FALSE)
  }
  mismatch <- names(expected)[!vapply(names(expected), function(name) {
    identical(artifact[[name]], expected[[name]])
  }, logical(1L))]
  if (length(mismatch)) {
    return(list(
      status = "stale_or_incomplete", reusable = FALSE,
      mismatch = mismatch
    ))
  }
  list(status = "complete", reusable = TRUE, artifact = artifact)
}
