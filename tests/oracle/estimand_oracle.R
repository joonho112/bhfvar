# Frozen probability-scale estimand oracle.
#
# This file uses base R only and consumes primitive vectors/matrices. It does
# not accept prepared package objects and does not call implementation helpers.

.oracle_abort <- function(message) {
  stop(message, call. = FALSE)
}

.oracle_inv_logit <- function(x) {
  positive <- x >= 0
  answer <- numeric(length(x))
  answer[positive] <- 1 / (1 + exp(-x[positive]))
  exp_x <- exp(x[!positive])
  answer[!positive] <- exp_x / (1 + exp_x)
  answer
}

.oracle_align_named <- function(x, labels, name, nonnegative = FALSE,
                                positive = FALSE) {
  if (!is.numeric(x) || length(x) != length(labels) || is.null(names(x)) ||
      anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x)) ||
      !setequal(names(x), labels) || anyNA(x) || any(!is.finite(x)) ||
      (nonnegative && any(x < 0)) || (positive && any(x <= 0))) {
    .oracle_abort(paste0(name, " must be a finite exact-set named vector"))
  }
  unname(x[labels])
}

.oracle_align_matrix <- function(x, n_draw, labels, draw_id, name) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != n_draw ||
      ncol(x) != length(labels) || is.null(colnames(x)) ||
      anyDuplicated(colnames(x)) || !setequal(colnames(x), labels) ||
      is.null(rownames(x)) ||
      !identical(as.character(rownames(x)), draw_id) ||
      anyNA(x) || any(!is.finite(x))) {
    .oracle_abort(paste0(name, " must be a finite draw-by-level matrix"))
  }
  x[, labels, drop = FALSE]
}

.oracle_summary <- function(state_probability, population_share) {
  mean_probability <- sum(population_share * state_probability)
  between <- sum(population_share *
                   (state_probability - mean_probability)^2)
  within <- sum(population_share * state_probability *
                  (1 - state_probability))
  total <- between + within
  proportion <- if (total > 0) between / total else 0
  c(
    mean = mean_probability,
    between = between,
    within = within,
    total = total,
    proportion = proportion
  )
}

bhf_reference_oracle <- function(design, draws, vhat = NULL,
                                 return_individual = FALSE) {
  required_design <- c(
    "domain_labels", "stratum_labels", "psu_labels", "state_id",
    "stratum_id", "psu_flat_id", "w_lik", "population_share"
  )
  if (!is.list(design) || any(!required_design %in% names(design))) {
    .oracle_abort("design is missing required primitive fields")
  }
  domain_labels <- as.character(design$domain_labels)
  stratum_labels <- as.character(design$stratum_labels)
  psu_labels <- as.character(design$psu_labels)
  if (!length(domain_labels) || !length(stratum_labels) || !length(psu_labels) ||
      anyNA(c(domain_labels, stratum_labels, psu_labels)) ||
      any(!nzchar(c(domain_labels, stratum_labels, psu_labels))) ||
      anyDuplicated(domain_labels) || anyDuplicated(stratum_labels) ||
      anyDuplicated(psu_labels)) {
    .oracle_abort("design labels must be nonblank and unique")
  }

  state_id <- design$state_id
  stratum_id <- design$stratum_id
  psu_flat_id <- design$psu_flat_id
  w_lik <- design$w_lik
  N <- length(state_id)
  S <- length(domain_labels)
  H <- length(stratum_labels)
  J <- length(psu_labels)
  valid_id <- function(x, upper) {
    is.numeric(x) && length(x) == N && !anyNA(x) && all(is.finite(x)) &&
      all(x == as.integer(x)) && all(x >= 1L & x <= upper) &&
      identical(sort(unique(as.integer(x))), seq_len(upper))
  }
  if (!valid_id(state_id, S) || !valid_id(stratum_id, H) ||
      !valid_id(psu_flat_id, J)) {
    .oracle_abort("design IDs must be complete consecutive observed levels")
  }
  state_id <- as.integer(state_id)
  stratum_id <- as.integer(stratum_id)
  psu_flat_id <- as.integer(psu_flat_id)
  psu_strata <- split(stratum_id, psu_flat_id)
  if (any(vapply(psu_strata, function(x) length(unique(x)) != 1L,
                 logical(1)))) {
    .oracle_abort("each psu_flat_id must be nested in exactly one stratum")
  }
  if (!is.numeric(w_lik) || length(w_lik) != N || anyNA(w_lik) ||
      any(!is.finite(w_lik)) || any(w_lik <= 0) ||
      abs(sum(w_lik) - N) > 1e-10 * max(1, N)) {
    .oracle_abort("w_lik must be positive, finite, and globally sum to N")
  }
  population_share <- .oracle_align_named(
    design$population_share, domain_labels, "population_share", positive = TRUE
  )
  if (abs(sum(population_share) - 1) > 1e-12) {
    .oracle_abort("population_share must sum to one")
  }
  vhat_aligned <- if (is.null(vhat)) NULL else .oracle_align_named(
    vhat, domain_labels, "vhat", nonnegative = TRUE
  )

  required_draws <- c(
    "alpha", "u_state", "u_stratum", "u_psu",
    "sigma_state", "sigma_stratum", "sigma_psu"
  )
  if (!is.list(draws) || any(!required_draws %in% names(draws)) ||
      !is.numeric(draws$alpha) || !length(draws$alpha) ||
      anyNA(draws$alpha) || any(!is.finite(draws$alpha))) {
    .oracle_abort("draws are missing finite alpha/effect/sigma fields")
  }
  D <- length(draws$alpha)
  draw_id <- if (is.null(draws$draw_id)) as.character(seq_len(D)) else {
    if (length(draws$draw_id) != D || anyNA(draws$draw_id) ||
        anyDuplicated(draws$draw_id)) {
      .oracle_abort("draw_id must be unique and aligned")
    }
    as.character(draws$draw_id)
  }
  u_state <- .oracle_align_matrix(
    draws$u_state, D, domain_labels, draw_id, "u_state"
  )
  u_stratum <- .oracle_align_matrix(
    draws$u_stratum, D, stratum_labels, draw_id, "u_stratum"
  )
  u_psu <- .oracle_align_matrix(
    draws$u_psu, D, psu_labels, draw_id, "u_psu"
  )
  sigmas <- lapply(
    c("sigma_state", "sigma_stratum", "sigma_psu"),
    function(name) {
      value <- draws[[name]]
      if (!is.numeric(value) || length(value) != D || anyNA(value) ||
          any(!is.finite(value)) || any(value < 0)) {
        .oracle_abort(paste0(name, " must be finite and nonnegative"))
      }
      value
    }
  )

  state_rows <- split(
    seq_len(N), factor(state_id, levels = seq_len(S))
  )
  state_weight <- vapply(
    state_rows, function(index) sum(w_lik[index]), numeric(1)
  )
  if (any(state_weight <= 0)) {
    .oracle_abort("every state must have positive likelihood-weight support")
  }

  summary_names <- c("mean", "between", "within", "total", "proportion")
  p_a <- matrix(NA_real_, D, S, dimnames = list(draw_id, domain_labels))
  p_b <- p_a
  a_summary <- matrix(NA_real_, D, 5L, dimnames = list(draw_id, summary_names))
  b_summary <- a_summary
  bernoulli_state <- p_a
  mixture_state <- p_a
  p_individual <- if (isTRUE(return_individual)) {
    row_labels <- if (is.null(design$row_id)) as.character(seq_len(N)) else
      as.character(design$row_id)
    if (length(row_labels) != N || anyNA(row_labels) || anyDuplicated(row_labels)) {
      .oracle_abort("row_id must be unique and aligned")
    }
    matrix(NA_real_, D, N, dimnames = list(draw_id, row_labels))
  } else NULL

  for (draw in seq_len(D)) {
    p_a[draw, ] <- .oracle_inv_logit(draws$alpha[draw] + u_state[draw, ])
    a_summary[draw, ] <- .oracle_summary(p_a[draw, ], population_share)

    eta <- draws$alpha[draw] + u_state[draw, state_id] +
      u_stratum[draw, stratum_id] + u_psu[draw, psu_flat_id]
    probability <- .oracle_inv_logit(eta)
    if (!is.null(p_individual)) p_individual[draw, ] <- probability

    for (state in seq_len(S)) {
      index <- state_rows[[state]]
      q <- w_lik[index] / state_weight[state]
      p_b[draw, state] <- sum(q * probability[index])
      bernoulli_state[draw, state] <-
        sum(q * probability[index] * (1 - probability[index]))
      mixture_state[draw, state] <-
        sum(q * (probability[index] - p_b[draw, state])^2)
    }
    mean_b <- sum(population_share * p_b[draw, ])
    between_b <- sum(population_share * (p_b[draw, ] - mean_b)^2)
    within_b <- sum(population_share *
                      (bernoulli_state[draw, ] + mixture_state[draw, ]))
    total_b <- between_b + within_b
    b_summary[draw, ] <- c(
      mean_b, between_b, within_b, total_b,
      if (total_b > 0) between_b / total_b else 0
    )
  }

  correction <- if (is.null(vhat_aligned)) NULL else
    sum(population_share * vhat_aligned)
  a_star <- list(
    available = FALSE,
    summary = NULL,
    correction = NULL,
    at_boundary = NULL,
    truncated = NULL,
    proportion_defined = NULL
  )
  if (!is.null(correction)) {
    a_star_summary <- a_summary
    corrected_between <- pmax(0, a_summary[, "between"] - correction)
    a_star_summary[, "between"] <- corrected_between
    a_star_summary[, "total"] <- corrected_between + a_summary[, "within"]
    a_star_summary[, "proportion"] <- ifelse(
      a_star_summary[, "total"] > 0,
      corrected_between / a_star_summary[, "total"],
      0
    )
    a_star <- list(
      available = TRUE,
      summary = a_star_summary,
      correction = rep(correction, D),
      at_boundary = corrected_between == 0,
      truncated = a_summary[, "between"] - correction < 0,
      proportion_defined = a_star_summary[, "total"] > 0
    )
  }

  variance_state <- sigmas[[1]]^2
  variance_stratum <- sigmas[[2]]^2
  variance_psu <- sigmas[[3]]^2
  variance_logistic <- rep(pi^2 / 3, D)
  latent_total <- variance_state + variance_stratum + variance_psu +
    variance_logistic
  latent <- cbind(
    var_state = variance_state,
    var_stratum = variance_stratum,
    var_psu = variance_psu,
    var_logistic = variance_logistic,
    total = latent_total,
    icc_state = variance_state / latent_total,
    icc_stratum = variance_stratum / latent_total,
    icc_psu = variance_psu / latent_total,
    icc_logistic = variance_logistic / latent_total
  )
  rownames(latent) <- draw_id

  gaps <- list(
    B_minus_A = b_summary - a_summary,
    A_minus_A_star = if (!isTRUE(a_star$available)) NULL else
      a_summary - a_star$summary
  )
  list(
    schema_version = "1.0.0",
    formula_authority = "math14030512-probability-decomposition-v1",
    assumptions = list(
      effects = "realized effects already satisfy article centering",
      psu = "psu_flat_id is nested in stratum_id",
      weights = "w_lik is globally mean-one and used once",
      universe = "every target state is observed",
      labels = "all named inputs use canonical exact-set alignment",
      draw_vectors = "alpha and sigma vectors are positional in draw_id order",
      total_zero = "numeric proportion is zero with defined flag false"
    ),
    dimensions = c(draws = D, observations = N, states = S,
                   strata = H, psus = J),
    draw_id = draw_id,
    latent = latent,
    A = list(
      p_state = p_a,
      summary = a_summary,
      proportion_defined = a_summary[, "total"] > 0
    ),
    A_star = a_star,
    B = list(
      p_state = p_b,
      within_bernoulli_state = bernoulli_state,
      within_mixture_state = mixture_state,
      within_bernoulli = rowSums(sweep(
        bernoulli_state, 2, population_share, `*`
      )),
      within_mixture = rowSums(sweep(
        mixture_state, 2, population_share, `*`
      )),
      summary = b_summary,
      proportion_defined = b_summary[, "total"] > 0,
      p_individual = p_individual
    ),
    gaps = gaps
  )
}
