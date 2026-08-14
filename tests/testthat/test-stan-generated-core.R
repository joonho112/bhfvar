# Static and deterministic contracts for Phase 5 generated quantities.
# Native Stan execution is deliberately reserved for the opt-in parity tier.

.bhf_gq_stan_path <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "inst", "stan", "bhf_hybrid.stan"),
    file.path("inst", "stan", "bhf_hybrid.stan"),
    system.file("stan", "bhf_hybrid.stan", package = "bhfvar")
  )
  existing <- unique(candidates[nzchar(candidates) & file.exists(candidates)])
  if (!length(existing)) stop("Cannot locate Stan source", call. = FALSE)
  normalizePath(existing[[1L]], mustWork = TRUE)
}

.bhf_gq_source <- function() {
  lines <- readLines(.bhf_gq_stan_path(), warn = FALSE)
  paste(sub("//.*$", "", lines), collapse = "\n")
}

.bhf_gq_block <- function(name) {
  source <- .bhf_gq_source()
  heading <- gsub(" ", "[[:space:]]+", name, fixed = TRUE)
  hit <- regexpr(
    paste0("(?m)^[[:space:]]*", heading, "[[:space:]]*\\{"),
    source,
    perl = TRUE
  )
  if (hit[[1L]] < 0L) stop("Cannot locate Stan block: ", name, call. = FALSE)
  matched <- regmatches(source, hit)
  open <- hit[[1L]] + regexpr("\\{", matched, perl = TRUE)[[1L]] - 1L
  chars <- strsplit(substring(source, open), "", fixed = TRUE)[[1L]]
  depth <- cumsum((chars == "{") - (chars == "}"))
  close <- which(depth == 0L)[[1L]]
  substring(source, open + 1L, open + close - 2L)
}

.bhf_gq_compact <- function(x) trimws(gsub("[[:space:]]+", " ", x))

.bhf_gq_manifest <- function() {
  lines <- strsplit(.bhf_gq_block("generated quantities"), "\n",
                    fixed = TRUE)[[1L]]
  declaration <- paste0(
    "^(?:real|int(?:<[^>]+>)?|vector(?:<[^>]+>)?\\[[^]]+\\]|",
    "array\\[[^]]+\\][[:space:]]+int(?:<[^>]+>)?)",
    "[[:space:]]+([A-Za-z][A-Za-z0-9_]*)[[:space:]]*;$"
  )
  answer <- character()
  depth <- 0L
  for (line in lines) {
    trimmed <- trimws(line)
    if (depth == 0L) {
      match <- regexec(declaration, trimmed, perl = TRUE)
      pieces <- regmatches(trimmed, match)[[1L]]
      if (length(pieces) == 2L) answer <- c(answer, pieces[[2L]])
    }
    chars <- strsplit(line, "", fixed = TRUE)[[1L]]
    depth <- depth + sum(chars == "{") - sum(chars == "}")
  }
  answer
}

.bhf_gq_region <- function(source, from, to) {
  start <- regexpr(from, source, fixed = TRUE)[[1L]]
  if (start < 0L) stop("Cannot locate region start: ", from, call. = FALSE)
  remainder <- substring(source, start)
  finish <- regexpr(to, remainder, fixed = TRUE)[[1L]]
  if (finish < 0L) stop("Cannot locate region end: ", to, call. = FALSE)
  substring(remainder, 1L, finish - 1L)
}

.bhf_gq_loop_table <- function(source) {
  matches <- gregexpr(
    "for[[:space:]]*\\([^)]+\\)[[:space:]]*\\{",
    source,
    perl = TRUE
  )[[1L]]
  if (identical(matches, -1L)) {
    return(data.frame(header = character(), depth = integer()))
  }
  lengths <- attr(matches, "match.length")
  headers <- substring(source, matches, matches + lengths - 1L)
  depths <- vapply(matches, function(position) {
    if (position == 1L) return(0L)
    chars <- strsplit(substring(source, 1L, position - 1L), "",
                      fixed = TRUE)[[1L]]
    sum(chars == "{") - sum(chars == "}")
  }, integer(1))
  data.frame(header = headers, depth = depths, stringsAsFactors = FALSE)
}

.bhf_gq_summary <- function(probability, share) {
  mean <- sum(share * probability)
  between <- sum(share * (probability - mean)^2)
  within <- sum(share * probability * (1 - probability))
  total <- between + within
  c(
    mean = mean,
    between = between,
    within = within,
    total = total,
    proportion = if (total > 0) between / total else 0,
    proportion_defined = as.numeric(total > 0)
  )
}

.bhf_gq_astar <- function(between, within, share, vhat,
                          use_deattenuation) {
  if (use_deattenuation) {
    correction <- sum(share * vhat)
    raw <- between - correction
    corrected <- max(0, raw)
    boundary <- corrected == 0
    truncated <- raw < 0
  } else {
    correction <- 0
    corrected <- between
    boundary <- FALSE
    truncated <- FALSE
  }
  total <- corrected + within
  c(
    correction = correction,
    between = corrected,
    within = within,
    total = total,
    proportion = if (total > 0) corrected / total else 0,
    at_boundary = as.numeric(boundary),
    truncated = as.numeric(truncated),
    proportion_defined = as.numeric(total > 0)
  )
}

.bhf_gq_b <- function(eta, state_id, w_lik, share) {
  probability <- stats::plogis(eta)
  S <- length(share)
  state_rows <- split(seq_along(state_id), factor(state_id, levels = seq_len(S)))
  state_weight <- vapply(state_rows, function(i) sum(w_lik[i]), numeric(1))
  p_state <- vapply(state_rows, function(i) {
    sum(w_lik[i] * probability[i]) / sum(w_lik[i])
  }, numeric(1))
  binomial_state <- vapply(state_rows, function(i) {
    sum(w_lik[i] * probability[i] * (1 - probability[i])) / sum(w_lik[i])
  }, numeric(1))
  mixture_state <- vapply(seq_len(S), function(s) {
    i <- state_rows[[s]]
    sum(w_lik[i] * (probability[i] - p_state[s])^2) / sum(w_lik[i])
  }, numeric(1))
  p_bar <- sum(share * p_state)
  within_binomial <- sum(share * binomial_state)
  within_mixture <- sum(share * mixture_state)
  within <- within_binomial + within_mixture
  between <- sum(share * (p_state - p_bar)^2)
  total <- between + within
  list(
    probability = probability,
    state_weight = state_weight,
    p_state = p_state,
    binomial_state = binomial_state,
    mixture_state = mixture_state,
    summary = c(
      mean = p_bar, between = between, within = within, total = total,
      proportion = if (total > 0) between / total else 0,
      proportion_defined = as.numeric(total > 0)
    ),
    within_binomial = within_binomial,
    within_mixture = within_mixture
  )
}

test_that("generated quantities expose exactly the approved output manifest", {
  expected <- c(
    "state_centering_residual", "stratum_centering_residual",
    "psu_centering_residual", "var_state_latent", "var_stratum_latent",
    "var_psu_latent", "var_level1_latent", "var_total_latent",
    "icc_state_latent", "icc_stratum_latent", "icc_psu_latent",
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
  expect_identical(.bhf_gq_manifest(), expected)
  expect_identical(length(expected), 55L)
})

test_that("latent diagnostics use three design variances and logistic level one", {
  gq <- .bhf_gq_compact(.bhf_gq_block("generated quantities"))
  patterns <- c(
    "var_state_latent = square\\(sigma_state\\)",
    "var_stratum_latent = square\\(sigma_stratum\\)",
    "var_psu_latent = square\\(sigma_psu\\)",
    "var_level1_latent = square\\(pi\\(\\)\\) / 3\\.0",
    paste0(
      "var_total_latent = var_state_latent \\+ var_stratum_latent \\+ ",
      "var_psu_latent \\+ var_level1_latent"
    ),
    "icc_state_latent = var_state_latent / var_total_latent",
    "icc_stratum_latent = var_stratum_latent / var_total_latent",
    "icc_psu_latent = var_psu_latent / var_total_latent"
  )
  for (pattern in patterns) expect_match(gq, pattern, perl = TRUE, info = pattern)

  variance <- c(state = 0.7^2, stratum = 0.4^2, psu = 0.25^2,
                level1 = pi^2 / 3)
  total <- sum(variance)
  icc <- variance[c("state", "stratum", "psu")] / total
  expect_true(all(is.finite(c(variance, total, icc))))
  expect_equal(sum(variance / total), 1, tolerance = 1e-15)
  expect_true(all(icc >= 0 & icc <= 1))
})

test_that("Estimand A is state-only and obeys weighted decomposition identities", {
  gq <- .bhf_gq_compact(.bhf_gq_block("generated quantities"))
  patterns <- c(
    "p_state_A\\[s\\] = inv_logit\\(alpha \\+ u_state\\[s\\]\\)",
    "p_bar_A = dot_product\\(w_state_pop_share, p_state_A\\)",
    paste0(
      "var_between_A \\+= w_state_pop_share\\[s\\] \\* ",
      "square\\(p_state_A\\[s\\] - p_bar_A\\)"
    ),
    paste0(
      "var_within_A \\+= w_state_pop_share\\[s\\] \\* p_state_A\\[s\\] \\* ",
      "\\(1\\.0 - p_state_A\\[s\\]\\)"
    ),
    "var_total_A = var_between_A \\+ var_within_A",
    "A_proportion_defined = var_total_A > 0",
    "prop_between_A = A_proportion_defined == 1"
  )
  for (pattern in patterns) expect_match(gq, pattern, perl = TRUE, info = pattern)
  expect_false(grepl(
    "p_state_A\\[s\\].*(u_stratum|u_psu|sigma_stratum|sigma_psu)",
    gq,
    perl = TRUE
  ))

  share <- c(0.5, 0.3, 0.2)
  p <- stats::plogis(-0.4 + c(0.3, -0.5, 0.1))
  a <- .bhf_gq_summary(p, share)
  expect_true(all(p > 0 & p < 1))
  expect_equal(a[["total"]], a[["between"]] + a[["within"]],
               tolerance = 1e-15)
  expect_equal(a[["total"]], a[["mean"]] * (1 - a[["mean"]]),
               tolerance = 1e-15)
  expect_true(a[["proportion_defined"]] == 1)
  zero <- .bhf_gq_summary(c(0, 0), c(0.5, 0.5))
  expect_identical(unname(zero[c("total", "proportion",
                                 "proportion_defined")]), c(0, 0, 0))
})

test_that("A-star changes only between variance and has exact boundary policy", {
  gq <- .bhf_gq_compact(.bhf_gq_block("generated quantities"))
  enabled <- c(
    "if \\(use_deattenuation == 1\\)",
    "mean_vhat_state = dot_product\\(w_state_pop_share, vhat_state\\)",
    "raw_between_A_star = var_between_A - mean_vhat_state",
    "var_between_A_star = fmax\\(0\\.0, raw_between_A_star\\)",
    "A_star_at_boundary = var_between_A_star == 0\\.0",
    "A_star_truncated = raw_between_A_star < 0\\.0"
  )
  disabled <- c(
    "mean_vhat_state = 0\\.0",
    "var_between_A_star = var_between_A",
    "A_star_at_boundary = 0", "A_star_truncated = 0"
  )
  identities <- c(
    "p_bar_A_star = p_bar_A",
    "var_within_A_star = var_within_A",
    "var_total_A_star = var_between_A_star \\+ var_within_A_star",
    "A_star_proportion_defined = var_total_A_star > 0",
    "prop_between_A_star = A_star_proportion_defined == 1"
  )
  for (pattern in c(enabled, disabled, identities)) {
    expect_match(gq, pattern, perl = TRUE, info = pattern)
  }
  expect_false(grepl("fmax\\(0\\.001", gq, perl = TRUE))

  share <- c(0.5, 0.5)
  disabled_case <- .bhf_gq_astar(0.25, 0.5, share, c(9, 9), FALSE)
  equal_case <- .bhf_gq_astar(0.25, 0.5, share, c(0.25, 0.25), TRUE)
  overshoot <- .bhf_gq_astar(0.25, 0.5, share, c(0.30, 0.30), TRUE)
  expect_equal(disabled_case[c("between", "within", "total", "proportion")],
               c(between = 0.25, within = 0.5, total = 0.75,
                 proportion = 1 / 3), tolerance = 1e-15)
  expect_identical(unname(disabled_case[c("correction", "at_boundary",
                                           "truncated")]), c(0, 0, 0))
  expect_identical(unname(equal_case[c("between", "at_boundary",
                                       "truncated")]), c(0, 1, 0))
  expect_identical(unname(overshoot[c("between", "at_boundary",
                                      "truncated")]), c(0, 1, 1))
})

test_that("Estimand B uses full eta and sequential O(N plus S) accumulators", {
  gq_raw <- .bhf_gq_block("generated quantities")
  gq <- .bhf_gq_compact(gq_raw)
  transformed <- .bhf_gq_compact(.bhf_gq_block("transformed data"))
  patterns <- c(
    "p = inv_logit\\(eta\\[i\\]\\)",
    "p_individual_B\\[i\\] = p",
    "state_weight_B\\[s\\] \\+= w_lik\\[i\\]",
    "p_state_B\\[s\\] \\+= w_lik\\[i\\] \\* p",
    "p_state_B\\[s\\] /= state_weight_B\\[s\\]",
    paste0(
      "within_binomial_state_B\\[s\\] \\+= w_lik\\[i\\] \\* p \\* ",
      "\\(1\\.0 - p\\)"
    ),
    paste0(
      "within_mixture_state_B\\[s\\] \\+= w_lik\\[i\\] \\* ",
      "square\\(p_individual_B\\[i\\] - p_state_B\\[s\\]\\)"
    )
  )
  for (pattern in patterns) expect_match(gq, pattern, perl = TRUE, info = pattern)
  expect_match(
    transformed,
    "state_weight_check\\[state_id\\[i\\]\\] \\+= w_lik\\[i\\]",
    perl = TRUE
  )
  expect_match(transformed, "min\\(state_weight_check\\) <= 0", perl = TRUE)

  region <- .bhf_gq_region(
    gq_raw,
    "state_weight_B = rep_vector(0.0, S);",
    "gap_B_minus_A_mean"
  )
  loops <- .bhf_gq_loop_table(region)
  expect_true(nrow(loops) > 0L)
  expect_true(all(loops$depth == 0L))
  expect_identical(sum(grepl("i in 1:N", loops$header, fixed = TRUE)), 2L)
  expect_identical(sum(grepl("s in 1:S", loops$header, fixed = TRUE)), 3L)
})

test_that("B binomial, mixture, total, proportions, and signed gaps close", {
  gq <- .bhf_gq_compact(.bhf_gq_block("generated quantities"))
  patterns <- c(
    paste0(
      "var_within_binomial_B = dot_product\\([[:space:]]*",
      "w_state_pop_share, within_binomial_state_B[[:space:]]*\\)"
    ),
    paste0(
      "var_within_mixture_B = dot_product\\([[:space:]]*",
      "w_state_pop_share, within_mixture_state_B[[:space:]]*\\)"
    ),
    "var_within_B = var_within_binomial_B \\+ var_within_mixture_B",
    "p_bar_B = dot_product\\(w_state_pop_share, p_state_B\\)",
    paste0(
      "var_between_B \\+= w_state_pop_share\\[s\\] \\* ",
      "square\\(p_state_B\\[s\\] - p_bar_B\\)"
    ),
    "var_total_B = var_between_B \\+ var_within_B",
    "B_proportion_defined = var_total_B > 0",
    "prop_between_B = B_proportion_defined == 1",
    "gap_B_minus_A_mean = p_bar_B - p_bar_A",
    "gap_B_minus_A_between = var_between_B - var_between_A",
    "gap_B_minus_A_within = var_within_B - var_within_A",
    "gap_B_minus_A_total = var_total_B - var_total_A",
    "gap_B_minus_A_proportion = prop_between_B - prop_between_A",
    "gap_A_minus_A_star_between = var_between_A - var_between_A_star",
    "gap_A_minus_A_star_total = var_total_A - var_total_A_star",
    "prop_between_A - prop_between_A_star"
  )
  for (pattern in patterns) expect_match(gq, pattern, perl = TRUE, info = pattern)

  state_id <- c(1L, 1L, 1L, 2L, 2L, 2L)
  w_lik <- c(0.4, 0.8, 1.2, 0.6, 1.0, 2.0)
  share <- c(0.65, 0.35)
  eta <- c(-2.0, -0.3, 1.2, -1.0, 0.6, 2.1)
  b <- .bhf_gq_b(eta, state_id, w_lik, share)
  expect_equal(sum(w_lik), length(w_lik), tolerance = 1e-15)
  expect_true(all(is.finite(unlist(b))))
  expect_true(all(b$probability > 0 & b$probability < 1))
  expect_true(b$within_binomial >= 0 && b$within_mixture >= 0)
  expect_equal(b$summary[["within"]],
               b$within_binomial + b$within_mixture, tolerance = 1e-15)
  expect_equal(b$summary[["total"]],
               b$summary[["between"]] + b$summary[["within"]],
               tolerance = 1e-15)
  expect_equal(b$summary[["total"]],
               b$summary[["mean"]] * (1 - b$summary[["mean"]]),
               tolerance = 1e-15)

  a <- .bhf_gq_summary(stats::plogis(c(-0.4, 0.2)), share)
  gap <- b$summary[names(a)] - a
  expect_equal(gap[["total"]], b$summary[["total"]] - a[["total"]],
               tolerance = 1e-15)
  expect_true(abs(gap[["mean"]]) > 1e-6)
})

test_that("finite state SD diagnostics are explicit, finite, and zero-safe", {
  gq <- .bhf_gq_compact(.bhf_gq_block("generated quantities"))
  patterns <- c(
    "sd_state_logit_unweighted = S > 1 \\? sd\\(u_state\\) : 0\\.0",
    paste0(
      "sd_state_probability_A_unweighted = S > 1 \\? ",
      "sd\\(p_state_A\\) : 0\\.0"
    ),
    "mean_logit_weighted = dot_product\\(w_state_pop_share, u_state\\)",
    paste0(
      "var_logit_weighted \\+= w_state_pop_share\\[s\\] \\* ",
      "square\\(u_state\\[s\\] - mean_logit_weighted\\)"
    ),
    paste0(
      "var_probability_A_weighted \\+= w_state_pop_share\\[s\\] \\* ",
      "square\\(p_state_A\\[s\\] - p_bar_A\\)"
    ),
    "sd_state_logit_weighted = sqrt\\(fmax\\(0\\.0, var_logit_weighted\\)\\)",
    paste0(
      "sd_state_probability_A_weighted = ",
      "sqrt\\(fmax\\(0\\.0, var_probability_A_weighted\\)\\)"
    )
  )
  for (pattern in patterns) expect_match(gq, pattern, perl = TRUE, info = pattern)

  u <- c(-0.8, 0.1, 0.7)
  p <- stats::plogis(-0.2 + u)
  share <- c(0.2, 0.5, 0.3)
  values <- c(
    logit_unweighted = stats::sd(u),
    probability_unweighted = stats::sd(p),
    logit_weighted = sqrt(sum(share * (u - sum(share * u))^2)),
    probability_weighted = sqrt(sum(share * (p - sum(share * p))^2)),
    singleton_logit = if (length(0.4) > 1) stats::sd(0.4) else 0,
    singleton_probability = if (length(0.7) > 1) stats::sd(0.7) else 0
  )
  expect_true(all(is.finite(values)))
  expect_true(all(values >= 0))
  expect_identical(unname(values[c("singleton_logit",
                                   "singleton_probability")]), c(0, 0))
})

test_that("obsolete Zeger, reliability, and ambiguous log-likelihood names stay absent", {
  source <- .bhf_gq_compact(.bhf_gq_source())
  forbidden <- c(
    "\\bc_factor\\b", "0\\.346", "\\breliability(?:_state|_avg)?\\b",
    "\\bp_state_marginal\\b", "\\bicc_prob\\b", "\\bicc_deatten\\b",
    "\\blog_lik\\s*\\["
  )
  for (pattern in forbidden) {
    expect_false(grepl(pattern, source, perl = TRUE), info = pattern)
  }
})

test_that("raw and pseudo log likelihoods and unweighted y-rep are unambiguous", {
  gq <- .bhf_gq_compact(.bhf_gq_block("generated quantities"))
  patterns <- c(
    "vector\\[N\\] log_lik_pseudo;",
    "vector\\[N\\] log_lik_raw;",
    "array\\[N\\] int<lower=0, upper=1> y_rep;",
    "log_lik_raw\\[i\\] = bernoulli_logit_lpmf\\(y\\[i\\] \\| eta\\[i\\]\\)",
    "log_lik_pseudo\\[i\\] = w_lik\\[i\\] \\* log_lik_raw\\[i\\]",
    "y_rep\\[i\\] = bernoulli_logit_rng\\(eta\\[i\\]\\)"
  )
  for (pattern in patterns) expect_match(gq, pattern, perl = TRUE, info = pattern)
  expect_false(grepl("y_rep\\[i\\][^;]*w_lik", gq, perl = TRUE))

  y <- c(0L, 1L, 1L, 0L)
  eta <- c(-1.2, -0.1, 0.8, 1.5)
  w_lik <- c(0.4, 0.8, 1.2, 1.6)
  raw <- stats::dbinom(y, 1, stats::plogis(eta), log = TRUE)
  pseudo <- w_lik * raw
  expect_true(all(is.finite(c(raw, pseudo))))
  expect_equal(pseudo, w_lik * raw, tolerance = 1e-15)
  expect_false(isTRUE(all.equal(raw, pseudo, tolerance = 1e-15)))
})
