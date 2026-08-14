# Deterministic article-aligned synthetic DGP used by validation fixtures.

.bhf_synthetic_abort <- function(message) {
  stop(message, call. = FALSE)
}

.bhf_synthetic_integer <- function(x, name, lower = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != as.integer(x) || x < lower || x > .Machine$integer.max) {
    .bhf_synthetic_abort(paste0(name, " must be one integer >= ", lower, "."))
  }
  as.integer(x)
}

.bhf_synthetic_with_seed <- function(seed, callback) {
  seed <- .bhf_synthetic_integer(seed, "seed", 1L)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  callback()
}

.bhf_synthetic_profile <- function(profile) {
  if (!is.character(profile) || length(profile) != 1L || is.na(profile) ||
      !profile %in% c("low", "high")) {
    .bhf_synthetic_abort("profile must be exactly 'low' or 'high'.")
  }
  if (identical(profile, "low")) {
    c(state = 0.3, stratum = 0.4, psu = 0.5)
  } else {
    c(state = 0.6, stratum = 0.7, psu = 0.9)
  }
}

.bhf_synthetic_center_state <- function(z, shares, sigma) {
  centered <- z - sum(shares * z)
  centered[[1L]] <- centered[[1L]] - sum(shares * centered) / shares[[1L]]
  sigma * centered
}

.bhf_synthetic_center_simple <- function(z, sigma) {
  centered <- z - mean(z)
  centered[[1L]] <- centered[[1L]] - sum(centered)
  sigma * centered
}

.bhf_synthetic_adler32 <- function(text) {
  bytes <- as.integer(charToRaw(enc2utf8(paste(text, collapse = "\n"))))
  a <- 1
  b <- 0
  for (byte in bytes) {
    a <- (a + byte) %% 65521
    b <- (b + a) %% 65521
  }
  paste0(sprintf("%04x", as.integer(b)), sprintf("%04x", as.integer(a)))
}

.bhf_synthetic_checksum <- function(x) {
  text <- utils::capture.output(base::dput(x))
  .bhf_synthetic_adler32(text)
}

.bhf_article_summary <- function(probability, shares) {
  mean_probability <- sum(shares * probability)
  between <- sum(shares * (probability - mean_probability)^2)
  within <- sum(shares * probability * (1 - probability))
  total <- between + within
  c(
    mean = mean_probability, between = between, within = within,
    total = total, proportion = if (total > 0) between / total else 0
  )
}

.bhf_article_truth_from_primitives <- function(data,
                                               alpha,
                                               u_state,
                                               u_stratum,
                                               u_psu,
                                               sigmas,
                                               population_shares,
                                               vhat_state) {
  required <- c("state_id", "stratum_id", "psu_flat_id", "w_lik")
  if (!is.data.frame(data) || !all(required %in% names(data)) || !nrow(data)) {
    .bhf_synthetic_abort("data is missing an article-truth primitive field.")
  }
  S <- length(u_state)
  H <- length(u_stratum)
  J <- length(u_psu)
  if (length(population_shares) != S || length(vhat_state) != S ||
      anyNA(c(alpha, u_state, u_stratum, u_psu, sigmas,
              population_shares, vhat_state)) ||
      any(!is.finite(c(alpha, u_state, u_stratum, u_psu, sigmas,
                       population_shares, vhat_state))) ||
      any(population_shares <= 0) ||
      abs(sum(population_shares) - 1) > 1e-12 || any(vhat_state < 0)) {
    .bhf_synthetic_abort("truth primitives are incomplete or invalid.")
  }
  state_id <- as.integer(data$state_id)
  stratum_id <- as.integer(data$stratum_id)
  psu_id <- as.integer(data$psu_flat_id)
  w_lik <- as.numeric(data$w_lik)
  if (any(state_id < 1L | state_id > S) ||
      any(stratum_id < 1L | stratum_id > H) ||
      any(psu_id < 1L | psu_id > J) || any(w_lik <= 0) ||
      abs(sum(w_lik) - nrow(data)) > 1e-10) {
    .bhf_synthetic_abort("truth data IDs or likelihood weights are invalid.")
  }

  p_state_A <- stats::plogis(alpha + u_state)
  A <- .bhf_article_summary(p_state_A, population_shares)
  correction <- sum(population_shares * vhat_state)
  raw_between_star <- unname(A[["between"]] - correction)
  between_star <- max(0, raw_between_star)
  total_star <- between_star + unname(A[["within"]])
  A_star <- c(
    mean = unname(A[["mean"]]), between = between_star,
    within = unname(A[["within"]]), total = total_star,
    proportion = if (total_star > 0) between_star / total_star else 0
  )

  eta <- alpha + u_state[state_id] + u_stratum[stratum_id] + u_psu[psu_id]
  p_individual_B <- stats::plogis(eta)
  state_weight <- numeric(S)
  p_state_B <- numeric(S)
  binomial_state <- numeric(S)
  for (s in seq_len(S)) {
    rows <- state_id == s
    state_weight[[s]] <- sum(w_lik[rows])
    p_state_B[[s]] <- sum(w_lik[rows] * p_individual_B[rows]) /
      state_weight[[s]]
    binomial_state[[s]] <- sum(
      w_lik[rows] * p_individual_B[rows] * (1 - p_individual_B[rows])
    ) / state_weight[[s]]
  }
  mixture_state <- vapply(seq_len(S), function(s) {
    rows <- state_id == s
    sum(w_lik[rows] * (p_individual_B[rows] - p_state_B[[s]])^2) /
      state_weight[[s]]
  }, numeric(1))
  p_bar_B <- sum(population_shares * p_state_B)
  between_B <- sum(population_shares * (p_state_B - p_bar_B)^2)
  within_binomial_B <- sum(population_shares * binomial_state)
  within_mixture_B <- sum(population_shares * mixture_state)
  within_B <- within_binomial_B + within_mixture_B
  total_B <- between_B + within_B
  B <- c(
    mean = p_bar_B, between = between_B,
    within_binomial = within_binomial_B,
    within_mixture = within_mixture_B, within = within_B,
    total = total_B,
    proportion = if (total_B > 0) between_B / total_B else 0
  )

  latent_variance <- c(
    state = sigmas[["state"]]^2,
    stratum = sigmas[["stratum"]]^2,
    psu = sigmas[["psu"]]^2,
    level1 = pi^2 / 3
  )
  latent_total <- sum(latent_variance)
  list(
    A = list(p_state = p_state_A, summary = A),
    A_star = list(
      correction = correction, summary = A_star,
      at_boundary = between_star == 0,
      truncated = raw_between_star < 0
    ),
    B = list(
      p_state = p_state_B, state_weight = state_weight,
      within_binomial_state = binomial_state,
      within_mixture_state = mixture_state, summary = B
    ),
    gaps = list(
      B_minus_A = B[c("mean", "between", "within", "total", "proportion")] -
        A[c("mean", "between", "within", "total", "proportion")],
      A_minus_A_star = A[c("between", "total", "proportion")] -
        A_star[c("between", "total", "proportion")]
    ),
    latent = c(
      latent_variance, total = latent_total,
      icc_state = latent_variance[["state"]] / latent_total,
      icc_stratum = latent_variance[["stratum"]] / latent_total,
      icc_psu = latent_variance[["psu"]] / latent_total
    )
  )
}

#' Generate a Small Article-Aligned Synthetic Survey
#'
#' Internal deterministic generator for scientific validation. It implements
#' the published low/high variance profiles and informative-weight tilt without
#' running the article's large Monte Carlo study.
#'
#' @keywords internal
generate_article_synthetic <- function(profile = c("low", "high"),
                                       rho = c(0, 0.5),
                                       n_states = 4L,
                                       n_strata = 3L,
                                       psus_per_stratum = 2L,
                                       observations_per_cell = 2L,
                                       alpha = -1.5,
                                       population_shares = NULL,
                                       vhat_state = NULL,
                                       effect_seed = 7106L,
                                       outcome_seed = 7107L) {
  if (identical(profile, c("low", "high"))) profile <- "low"
  if (identical(rho, c(0, 0.5))) rho <- 0
  sigmas <- .bhf_synthetic_profile(profile)
  if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) ||
      !rho %in% c(0, 0.5)) {
    .bhf_synthetic_abort("rho must be exactly 0 or 0.5.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha != -1.5) {
    .bhf_synthetic_abort("alpha must be exactly -1.5 for the article DGP.")
  }
  S <- .bhf_synthetic_integer(n_states, "n_states", 2L)
  H <- .bhf_synthetic_integer(n_strata, "n_strata", 2L)
  observations_per_cell <- .bhf_synthetic_integer(
    observations_per_cell, "observations_per_cell", 1L
  )
  if (length(psus_per_stratum) == 1L) {
    psus_per_stratum <- rep(psus_per_stratum, H)
  }
  if (!is.numeric(psus_per_stratum) || length(psus_per_stratum) != H ||
      anyNA(psus_per_stratum) || any(!is.finite(psus_per_stratum)) ||
      any(psus_per_stratum != as.integer(psus_per_stratum)) ||
      any(psus_per_stratum < 2L)) {
    .bhf_synthetic_abort(
      "psus_per_stratum must provide at least two integer PSUs per stratum."
    )
  }
  psus_per_stratum <- as.integer(psus_per_stratum)
  J <- sum(psus_per_stratum)
  state_labels <- sprintf("state%02d", seq_len(S))
  stratum_labels <- sprintf("stratum%02d", seq_len(H))
  if (is.null(population_shares)) population_shares <- seq_len(S)
  if (!is.numeric(population_shares) || length(population_shares) != S ||
      anyNA(population_shares) || any(!is.finite(population_shares)) ||
      any(population_shares <= 0)) {
    .bhf_synthetic_abort("population_shares must be S positive finite values.")
  }
  population_shares <- population_shares / sum(population_shares)
  names(population_shares) <- state_labels
  if (is.null(vhat_state)) vhat_state <- rep(0, S)
  if (!is.numeric(vhat_state) || length(vhat_state) != S ||
      anyNA(vhat_state) || any(!is.finite(vhat_state)) || any(vhat_state < 0)) {
    .bhf_synthetic_abort("vhat_state must be S nonnegative finite values.")
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
    block <- which(psu_stratum_id == h)
    u_psu[block] <- .bhf_synthetic_center_simple(
      z$psu[block], sigmas[["psu"]]
    )
  }
  names(u_state) <- state_labels
  names(u_stratum) <- stratum_labels
  psu_labels <- sprintf("psu%02d_%02d", psu_stratum_id, psu_local_id)
  names(u_psu) <- psu_labels

  data <- expand.grid(
    replicate_id = seq_len(observations_per_cell),
    state_id = seq_len(S), psu_flat_id = seq_len(J),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  data$stratum_id <- psu_stratum_id[data$psu_flat_id]
  data$psu_in_stratum_id <- psu_local_id[data$psu_flat_id]
  data$state <- state_labels[data$state_id]
  data$stratum <- stratum_labels[data$stratum_id]
  data$psu <- psu_labels[data$psu_flat_id]
  data$eta <- alpha + u_state[data$state_id] + u_stratum[data$stratum_id] +
    u_psu[data$psu_flat_id]
  data$probability <- stats::plogis(data$eta)
  data$raw_weight <- exp(-rho * u_stratum[data$stratum_id])
  data$weight <- data$raw_weight
  data$w_lik <- data$raw_weight * nrow(data) / sum(data$raw_weight)
  data$outcome <- .bhf_synthetic_with_seed(
    outcome_seed,
    function() stats::rbinom(nrow(data), size = 1L, prob = data$probability)
  )
  data$row_id <- seq_len(nrow(data))
  data <- data[c(
    "row_id", "outcome", "state", "stratum", "psu", "weight",
    "raw_weight", "w_lik", "state_id", "stratum_id",
    "psu_in_stratum_id", "psu_flat_id", "replicate_id", "eta",
    "probability"
  )]

  estimands <- .bhf_article_truth_from_primitives(
    data = data, alpha = alpha, u_state = u_state,
    u_stratum = u_stratum, u_psu = u_psu, sigmas = sigmas,
    population_shares = population_shares, vhat_state = vhat_state
  )
  truth <- list(
    schema_version = "article-synthetic-truth-1.0.0",
    authority = "Mathematics 14(3), 512 (2026), simulation design",
    config = list(
      profile = profile, rho = rho, alpha = alpha,
      effect_seed = as.integer(effect_seed), outcome_seed = as.integer(outcome_seed),
      dimensions = c(N = nrow(data), S = S, H = H, J = J),
      psus_per_stratum = psus_per_stratum,
      observations_per_cell = observations_per_cell
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
      }, numeric(1))
    ),
    weights = list(
      tilt = exp(-rho * u_stratum), raw_sum = sum(data$raw_weight),
      likelihood_sum = sum(data$w_lik), normalization = "global_mean_one"
    ),
    estimands = estimands
  )
  truth$checksums <- list(
    algorithm = "adler32-canonical-dput",
    data = .bhf_synthetic_checksum(data),
    truth_without_checksums = .bhf_synthetic_checksum(truth),
    outcomes = .bhf_synthetic_checksum(data$outcome)
  )
  structure(list(data = data, truth = truth),
            class = c("bhf_article_synthetic", "list"))
}
