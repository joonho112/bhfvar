# Prospective Gate G4-R validation helpers.

.bhf_g4r_abort <- function(message) {
  stop(structure(
    list(message = message, call = NULL),
    class = c("bhf_g4r_preflight_error", "error", "condition")
  ))
}

#' Frozen Gate G4-R Protocol
#' @noRd
bhf_g4r_protocol <- function() {
  protocol <- list(
    schema_version = "bhfvar-g4r-protocol-1.0.0",
    authority = list(
      protocol_document = "frozen simulation protocol (project records)",
      protocol_sha256 =
        "5c36a99403af3a7dbd3dfdaacd6a9a3d8249211f203ca2c161276ba91ce8df68",
      original_gate_status = "fail"
    ),
    dgp = list(
      rho = 0.5, alpha = -1.5, n_states = 20L, n_strata = 8L,
      psus_per_stratum = 3L, observations_per_cell = 2L,
      expected_n = 960L
    ),
    gap = list(
      equivalence_truth_max = 0.01,
      posterior_rope_max = 0.05,
      detectable_truth_min = 0.10,
      candidates_per_profile = 5000L,
      detectable_per_sign = 3L,
      equivalence_per_profile = 6L
    ),
    vhat = list(correction_fraction = 0.25, shape_min = 0.75,
                shape_max = 1.25),
    weight = list(kappa = 1.0, max_abs_delta_min = 0.10,
                  mean_abs_delta_min = 0.02),
    mcmc = list(
      chains = 4L, iter = 2000L, warmup = 1000L,
      kept_per_chain = 1000L, adapt_delta = 0.95,
      max_treedepth = 12L
    ),
    thresholds = list(
      divergences = 0L, treedepth_hits = 0L, ebfmi_min = 0.20,
      core_rhat_max = 1.01, ess_bulk_min = 400,
      ess_tail_min = 400, raw_effect_rhat_max = 1.05,
      interval_probability = 0.90, overall_coverage_min = 0.80,
      target_family_coverage_min = 20L,
      standardized_error_median_max = 1.0,
      standardized_error_q95_max = 2.5,
      regime_probability_min = 0.80,
      regime_pass_per_profile_min = 5L,
      detectable_pass_per_sign_min = 2L
    ),
    rerun = list(
      waves_max = 1L, reruns_per_job_max = 1L,
      eligibility = "sampler_diagnostic_failure_only",
      multiplier = 2L
    )
  )
  protocol$config_hash <- .bhf_synthetic_checksum(protocol)
  protocol
}

#' Calculate the Frozen G4-R Gap Statistic
#' @noRd
bhf_g4r_gap_metrics <- function(estimands,
                                protocol = bhf_g4r_protocol()) {
  if (!is.list(estimands) || is.null(estimands$A$summary) ||
      is.null(estimands$B$summary)) {
    .bhf_g4r_abort("estimands must contain A and B summary truth.")
  }
  a_between <- unname(estimands$A$summary[["between"]])
  b_between <- unname(estimands$B$summary[["between"]])
  if (length(a_between) != 1L || length(b_between) != 1L ||
      !is.finite(a_between) || !is.finite(b_between) || a_between < 0 ||
      b_between < 0) {
    .bhf_g4r_abort("A/B between truth must be finite nonnegative scalars.")
  }
  delta <- a_between - b_between
  ratio <- abs(delta) / max(a_between, 1e-12)
  gap <- protocol$gap
  regime <- if (ratio >= gap$detectable_truth_min) {
    "detectable"
  } else if (ratio <= gap$equivalence_truth_max) {
    "equivalence"
  } else {
    "grey"
  }
  direction <- if (delta > 0) "positive" else if (delta < 0) "negative" else
    "zero"
  list(
    delta_truth = delta, r_truth = ratio,
    regime = regime, direction = direction,
    a_between_truth = a_between, b_between_truth = b_between
  )
}

#' Map a Prospective Candidate Index to Frozen Seeds
#' @noRd
bhf_g4r_candidate_seeds <- function(profile, candidate_index) {
  if (!is.character(profile) || length(profile) != 1L ||
      !profile %in% c("low", "high")) {
    .bhf_g4r_abort("profile must be exactly 'low' or 'high'.")
  }
  limit <- bhf_g4r_protocol()$gap$candidates_per_profile
  if (!is.numeric(candidate_index) || length(candidate_index) != 1L ||
      is.na(candidate_index) || !is.finite(candidate_index) ||
      candidate_index != as.integer(candidate_index) || candidate_index < 1L ||
      candidate_index > limit) {
    .bhf_g4r_abort("candidate_index is outside the frozen candidate range.")
  }
  base <- if (identical(profile, "low")) 7310000L else 7410000L
  effect_seed <- base + 10L * as.integer(candidate_index)
  c(
    effect_seed = effect_seed,
    outcome_seed = effect_seed + 1L,
    mcmc_seed = effect_seed + 2L
  )
}

#' Fast Outcome-Free Gap Oracle for the Frozen Candidate Scan
#' @noRd
bhf_g4r_candidate_gap_oracle <- function(profile, effect_seed,
                                         protocol = bhf_g4r_protocol()) {
  dgp <- protocol$dgp
  sigmas <- .bhf_synthetic_profile(profile)
  effect_seed <- .bhf_synthetic_integer(effect_seed, "effect_seed", 1L)
  S <- dgp$n_states
  H <- dgp$n_strata
  psus_per_stratum <- rep(dgp$psus_per_stratum, H)
  J <- sum(psus_per_stratum)
  shares <- seq_len(S) / sum(seq_len(S))
  z <- .bhf_synthetic_with_seed(effect_seed, function() {
    list(
      state = stats::rnorm(S), stratum = stats::rnorm(H),
      psu = stats::rnorm(J)
    )
  })
  u_state <- .bhf_synthetic_center_state(
    z$state, shares, sigmas[["state"]]
  )
  u_stratum <- .bhf_synthetic_center_simple(
    z$stratum, sigmas[["stratum"]]
  )
  psu_stratum_id <- rep(seq_len(H), psus_per_stratum)
  u_psu <- numeric(J)
  for (h in seq_len(H)) {
    rows <- which(psu_stratum_id == h)
    u_psu[rows] <- .bhf_synthetic_center_simple(
      z$psu[rows], sigmas[["psu"]]
    )
  }
  p_A <- stats::plogis(dgp$alpha + u_state)
  p_bar_A <- sum(shares * p_A)
  a_between <- sum(shares * (p_A - p_bar_A)^2)
  design_effect <- u_stratum[psu_stratum_id] + u_psu
  probability <- stats::plogis(
    dgp$alpha + outer(u_state, design_effect, "+")
  )
  weight_psu <- exp(-dgp$rho * u_stratum[psu_stratum_id])
  p_B <- as.numeric(probability %*% weight_psu / sum(weight_psu))
  p_bar_B <- sum(shares * p_B)
  b_between <- sum(shares * (p_B - p_bar_B)^2)
  estimands <- list(
    A = list(summary = c(between = a_between)),
    B = list(summary = c(between = b_between))
  )
  c(
    bhf_g4r_gap_metrics(estimands, protocol),
    list(effect_seed = effect_seed)
  )
}

#' Select Frozen G4-R Regime Quotas from an Outcome-Blind Ledger
#' @noRd
bhf_g4r_select_candidates <- function(ledger,
                                      protocol = bhf_g4r_protocol()) {
  required <- c(
    "profile", "candidate_index", "effect_seed", "delta_truth", "r_truth",
    "regime", "direction"
  )
  forbidden <- intersect(
    names(ledger), c("outcome", "posterior", "draws", "fit", "stanfit")
  )
  if (!is.data.frame(ledger) || !all(required %in% names(ledger)) ||
      length(forbidden)) {
    .bhf_g4r_abort(
      "selector ledger is incomplete or contains outcome/posterior fields."
    )
  }
  if (anyDuplicated(ledger[c("profile", "candidate_index")]) ||
      any(!ledger$profile %in% c("low", "high")) ||
      anyNA(ledger[required])) {
    .bhf_g4r_abort("selector ledger keys and values must be complete and unique.")
  }
  ledger <- ledger[order(
    match(ledger$profile, c("low", "high")), ledger$candidate_index
  ), , drop = FALSE]
  quota_sign <- protocol$gap$detectable_per_sign
  quota_equiv <- protocol$gap$equivalence_per_profile
  selected <- list()
  position <- 0L
  for (profile in c("low", "high")) {
    block <- ledger$profile == profile
    groups <- list(
      detectable_positive = which(
        block & ledger$regime == "detectable" &
          ledger$direction == "positive"
      ),
      detectable_negative = which(
        block & ledger$regime == "detectable" &
          ledger$direction == "negative"
      ),
      equivalence = which(block & ledger$regime == "equivalence")
    )
    quotas <- c(quota_sign, quota_sign, quota_equiv)
    names(quotas) <- names(groups)
    for (group in names(groups)) {
      rows <- groups[[group]]
      if (length(rows) < quotas[[group]]) {
        .bhf_g4r_abort(paste0(
          "frozen candidate pool did not fill ", profile, "/", group,
          " quota."
        ))
      }
      keep <- rows[seq_len(quotas[[group]])]
      value <- ledger[keep, , drop = FALSE]
      value$selection_group <- group
      value$selection_rank <- seq_len(nrow(value))
      value$job_id <- sprintf(
        "g4r-%s-%s-%02d", profile,
        gsub("detectable_", "det-", gsub("equivalence", "equiv", group)),
        value$selection_rank
      )
      seeds <- t(vapply(seq_len(nrow(value)), function(index) {
        bhf_g4r_candidate_seeds(profile, value$candidate_index[[index]])
      }, numeric(3L)))
      value$outcome_seed <- as.integer(seeds[, "outcome_seed"])
      value$mcmc_seed <- as.integer(seeds[, "mcmc_seed"])
      position <- position + 1L
      selected[[position]] <- value
    }
  }
  result <- do.call(rbind, selected)
  rownames(result) <- NULL
  result
}

#' Select the Approved G4-R1 Sign-Agnostic Regime Quotas
#' @noRd
bhf_g4r1_select_candidates <- function(ledger,
                                       protocol = bhf_g4r_protocol()) {
  required <- c(
    "profile", "candidate_index", "effect_seed", "delta_truth", "r_truth",
    "regime", "direction"
  )
  forbidden <- intersect(
    names(ledger), c("outcome", "posterior", "draws", "fit", "stanfit")
  )
  if (!is.data.frame(ledger) || !all(required %in% names(ledger)) ||
      length(forbidden) || anyDuplicated(ledger[c("profile","candidate_index")]) ||
      anyNA(ledger[required]) || any(!ledger$profile %in% c("low","high"))) {
    .bhf_g4r_abort("G4-R1 ledger is incomplete, duplicated, or outcome-aware.")
  }
  ledger <- ledger[order(
    match(ledger$profile,c("low","high")), ledger$candidate_index
  ),,drop=FALSE]
  selected <- list()
  position <- 0L
  for (profile in c("low","high")) {
    block <- ledger[ledger$profile==profile,,drop=FALSE]
    groups <- list(
      detectable=block[block$regime=="detectable",,drop=FALSE],
      equivalence=block[block$regime=="equivalence",,drop=FALSE]
    )
    quotas <- c(
      detectable=2L*protocol$gap$detectable_per_sign,
      equivalence=protocol$gap$equivalence_per_profile
    )
    for (group in names(groups)) {
      value <- groups[[group]]
      if (nrow(value)<quotas[[group]]) {
        .bhf_g4r_abort(paste0(
          "G4-R1 candidate pool did not fill ",profile,"/",group," quota."
        ))
      }
      value <- value[seq_len(quotas[[group]]),,drop=FALSE]
      value$selection_group <- group
      value$selection_rank <- seq_len(nrow(value))
      value$job_id <- sprintf(
        "g4r1-%s-%s-%02d",profile,
        ifelse(group=="detectable","det","equiv"),value$selection_rank
      )
      seeds <- t(vapply(seq_len(nrow(value)),function(index) {
        bhf_g4r_candidate_seeds(profile,value$candidate_index[[index]])
      },numeric(3L)))
      value$outcome_seed <- as.integer(seeds[,"outcome_seed"])
      value$mcmc_seed <- as.integer(seeds[,"mcmc_seed"])
      position <- position+1L
      selected[[position]] <- value
    }
  }
  result <- do.call(rbind,selected)
  rownames(result) <- NULL
  result
}

#' Adjudicate One Frozen G4-R1 Gap Regime from Posterior Draws
#' @noRd
bhf_g4r1_regime_evidence <- function(delta_draw, a_between_draw,
                                     delta_truth, regime,
                                     protocol=bhf_g4r_protocol()) {
  if (!is.numeric(delta_draw) || !is.numeric(a_between_draw) ||
      length(delta_draw)<1L || length(delta_draw)!=length(a_between_draw) ||
      anyNA(delta_draw) || anyNA(a_between_draw) ||
      any(!is.finite(delta_draw)) || any(!is.finite(a_between_draw)) ||
      any(a_between_draw<0) || !is.numeric(delta_truth) ||
      length(delta_truth)!=1L || !is.finite(delta_truth) ||
      !is.character(regime) || length(regime)!=1L ||
      !regime %in% c("detectable","equivalence")) {
    .bhf_g4r_abort("regime evidence inputs are incomplete or invalid.")
  }
  threshold <- protocol$thresholds$regime_probability_min
  if (regime=="detectable") {
    if (delta_truth==0) {
      .bhf_g4r_abort("detectable truth must have a nonzero sign.")
    }
    probability <- mean(sign(delta_draw)==sign(delta_truth))
    metric <- "correct_sign_probability"
  } else {
    ratio <- abs(delta_draw)/pmax(a_between_draw,1e-12)
    probability <- mean(ratio<=protocol$gap$posterior_rope_max)
    metric <- "practical_equivalence_probability"
  }
  list(
    regime=regime, metric=metric, probability=probability,
    threshold=threshold, pass=probability>=threshold
  )
}

#' Generate an Outcome-Free Expanded G4-R Truth
#' @noRd
generate_g4r_truth_only <- function(profile, effect_seed,
                                    vhat_state = NULL) {
  protocol <- bhf_g4r_protocol()
  dgp <- protocol$dgp
  sigmas <- .bhf_synthetic_profile(profile)
  effect_seed <- .bhf_synthetic_integer(effect_seed, "effect_seed", 1L)
  S <- dgp$n_states
  H <- dgp$n_strata
  psus_per_stratum <- rep(dgp$psus_per_stratum, H)
  J <- sum(psus_per_stratum)
  state_labels <- sprintf("state%02d", seq_len(S))
  stratum_labels <- sprintf("stratum%02d", seq_len(H))
  population_shares <- seq_len(S) / sum(seq_len(S))
  names(population_shares) <- state_labels
  if (is.null(vhat_state)) vhat_state <- rep(0, S)
  if (!is.numeric(vhat_state) || length(vhat_state) != S ||
      anyNA(vhat_state) || any(!is.finite(vhat_state)) ||
      any(vhat_state < 0)) {
    .bhf_g4r_abort("vhat_state must be S finite nonnegative values.")
  }
  names(vhat_state) <- state_labels

  z <- .bhf_synthetic_with_seed(effect_seed, function() {
    list(
      state = stats::rnorm(S), stratum = stats::rnorm(H),
      psu = stats::rnorm(J)
    )
  })
  u_state <- .bhf_synthetic_center_state(
    z$state, population_shares, sigmas[["state"]]
  )
  u_stratum <- .bhf_synthetic_center_simple(
    z$stratum, sigmas[["stratum"]]
  )
  psu_stratum_id <- rep(seq_len(H), psus_per_stratum)
  psu_local_id <- sequence(psus_per_stratum)
  u_psu <- numeric(J)
  for (h in seq_len(H)) {
    rows <- which(psu_stratum_id == h)
    u_psu[rows] <- .bhf_synthetic_center_simple(
      z$psu[rows], sigmas[["psu"]]
    )
  }
  names(u_state) <- state_labels
  names(u_stratum) <- stratum_labels
  psu_labels <- sprintf("psu%02d_%02d", psu_stratum_id, psu_local_id)
  names(u_psu) <- psu_labels

  data <- expand.grid(
    replicate_id = seq_len(dgp$observations_per_cell),
    state_id = seq_len(S), psu_flat_id = seq_len(J),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  data$stratum_id <- psu_stratum_id[data$psu_flat_id]
  data$psu_in_stratum_id <- psu_local_id[data$psu_flat_id]
  data$state <- state_labels[data$state_id]
  data$stratum <- stratum_labels[data$stratum_id]
  data$psu <- psu_labels[data$psu_flat_id]
  data$eta <- dgp$alpha + u_state[data$state_id] +
    u_stratum[data$stratum_id] + u_psu[data$psu_flat_id]
  data$probability <- stats::plogis(data$eta)
  data$raw_weight <- exp(-dgp$rho * u_stratum[data$stratum_id])
  data$weight <- data$raw_weight
  data$w_lik <- data$raw_weight * nrow(data) / sum(data$raw_weight)
  data$row_id <- seq_len(nrow(data))
  data <- data[c(
    "row_id", "state", "stratum", "psu", "weight", "raw_weight",
    "w_lik", "state_id", "stratum_id", "psu_in_stratum_id",
    "psu_flat_id", "replicate_id", "eta", "probability"
  )]
  if (nrow(data) != dgp$expected_n) {
    .bhf_g4r_abort("expanded truth-only DGP did not produce frozen N=960.")
  }

  estimands <- .bhf_article_truth_from_primitives(
    data = data, alpha = dgp$alpha, u_state = u_state,
    u_stratum = u_stratum, u_psu = u_psu, sigmas = sigmas,
    population_shares = population_shares, vhat_state = vhat_state
  )
  truth <- list(
    schema_version = "bhfvar-g4r-truth-only-1.0.0",
    authority = "frozen prospective simulation protocol",
    outcome_generated = FALSE,
    config = list(
      profile = profile, rho = dgp$rho, alpha = dgp$alpha,
      effect_seed = effect_seed,
      dimensions = c(N = nrow(data), S = S, H = H, J = J),
      psus_per_stratum = psus_per_stratum,
      observations_per_cell = dgp$observations_per_cell
    ),
    sigmas = sigmas,
    population_shares = population_shares,
    vhat_state = vhat_state,
    effects = list(state = u_state, stratum = u_stratum, psu = u_psu),
    psu_structure = data.frame(
      psu = psu_labels, psu_flat_id = seq_len(J),
      stratum_id = psu_stratum_id,
      psu_in_stratum_id = psu_local_id,
      stringsAsFactors = FALSE
    ),
    centering_residuals = list(
      state_weighted = sum(population_shares * u_state),
      stratum_mean = mean(u_stratum),
      psu_within_stratum = vapply(seq_len(H), function(h) {
        mean(u_psu[psu_stratum_id == h])
      }, numeric(1L))
    ),
    weights = list(
      mode = "article_baseline",
      raw_sum = sum(data$raw_weight), likelihood_sum = sum(data$w_lik),
      normalization = "global_mean_one"
    ),
    estimands = estimands,
    gap = bhf_g4r_gap_metrics(estimands, protocol)
  )
  truth$checksums <- list(
    algorithm = "adler32-canonical-dput",
    data_without_outcome = .bhf_synthetic_checksum(data),
    truth_without_checksums = .bhf_synthetic_checksum(truth)
  )
  structure(
    list(data = data, truth = truth),
    class = c("bhf_g4r_truth_only", "list")
  )
}

#' Construct the Frozen Nonzero G4-R Sampling-Variance Truth
#' @noRd
bhf_g4r_vhat_from_truth <- function(truth,
                                    protocol = bhf_g4r_protocol()) {
  if (!is.list(truth) || is.null(truth$population_shares) ||
      is.null(truth$estimands$A$summary)) {
    .bhf_g4r_abort("truth is missing shares or A summary for vhat construction.")
  }
  shares <- truth$population_shares
  a_between <- unname(truth$estimands$A$summary[["between"]])
  if (!is.numeric(shares) || anyNA(shares) || any(shares <= 0) ||
      abs(sum(shares) - 1) > 1e-12 || !is.finite(a_between) ||
      a_between <= 0) {
    .bhf_g4r_abort("truth shares and positive A-between are required.")
  }
  S <- length(shares)
  shape <- seq(
    protocol$vhat$shape_min, protocol$vhat$shape_max, length.out = S
  )
  normalized_shape <- shape / sum(shares * shape)
  correction <- protocol$vhat$correction_fraction * a_between
  values <- correction * normalized_shape
  names(values) <- names(shares)
  list(
    values = values,
    target_correction = correction,
    observed_correction = sum(shares * values),
    correction_fraction = protocol$vhat$correction_fraction,
    normalized_shape = stats::setNames(normalized_shape, names(shares)),
    source = "supplied_simulation_truth"
  )
}

#' Generate an Outcome-Free G4-R Truth with Nonzero A-star Correction
#' @noRd
generate_g4r_truth_with_vhat <- function(profile, effect_seed) {
  zero <- generate_g4r_truth_only(profile, effect_seed)
  design <- bhf_g4r_vhat_from_truth(zero$truth)
  result <- generate_g4r_truth_only(profile, effect_seed, design$values)
  if (!identical(zero$data, result$data)) {
    .bhf_g4r_abort("adding vhat changed outcome-free DGP rows or weights.")
  }
  result$truth$vhat_design <- design
  result$truth$checksums$truth_without_checksums <- NULL
  result$truth$checksums$truth_without_checksums <- .bhf_synthetic_checksum(
    result$truth[names(result$truth) != "checksums"]
  )
  result
}

#' Materialize the Frozen Nondegenerate G4-R Weight Arm
#' @noRd
bhf_g4r_weight_arm <- function(truth_only,
                               protocol = bhf_g4r_protocol()) {
  if (!inherits(truth_only, "bhf_g4r_truth_only") ||
      is.null(truth_only$truth$effects$stratum)) {
    .bhf_g4r_abort("truth_only must be a complete bhf_g4r_truth_only object.")
  }
  data <- truth_only$data
  truth <- truth_only$truth
  S <- truth$config$dimensions[["S"]]
  H <- truth$config$dimensions[["H"]]
  g <- seq(-1, 1, length.out = S)
  h <- seq(-1, 1, length.out = H)
  kappa <- protocol$weight$kappa
  raw <- exp(
    -protocol$dgp$rho * truth$effects$stratum[data$stratum_id] +
      kappa * g[data$state_id] * h[data$stratum_id]
  )
  mean_one <- scale_likelihood_weights(raw, "mean_one", data$state_id)
  warning_messages <- character()
  legacy <- withCallingHandlers(
    scale_likelihood_weights(raw, "legacy_d2", data$state_id),
    warning = function(condition) {
      warning_messages <<- c(warning_messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  difference <- mean_one - legacy
  metrics <- list(
    max_abs_delta = max(abs(difference)),
    mean_abs_delta = mean(abs(difference)),
    mean_one_sum = sum(mean_one), legacy_d2_sum = sum(legacy),
    contrast_pass = max(abs(difference)) >=
      protocol$weight$max_abs_delta_min &&
      mean(abs(difference)) >= protocol$weight$mean_abs_delta_min
  )
  if (!metrics$contrast_pass) {
    .bhf_g4r_abort("frozen weight arm did not satisfy nondegenerate contrast.")
  }
  data$weight <- raw
  data$raw_weight <- raw
  data$w_lik <- mean_one
  list(
    schema_version = "bhfvar-g4r-weight-arm-1.0.0",
    data = data, raw_weight = raw, mean_one = mean_one,
    legacy_d2 = legacy, difference = difference, metrics = metrics,
    design = list(
      kappa = kappa, state_score = g, stratum_score = h,
      formula = "exp(-0.5*u_stratum[h] + kappa*g_s*h_h)"
    ),
    warning_messages = unique(warning_messages)
  )
}
