# Sensitivity of Estimand A* to the supplied sampling variances.
#
# A* subtracts a fixed, weighted sampling-variance correction from Estimand A's
# between-domain variance (article Equation 19). The V_s are treated as known,
# so their estimation uncertainty never reaches the posterior. This diagnostic
# does not change that; it shows how far A* would move if the V_s were off by a
# given factor, so the reader can see what the fixed-input assumption is worth.
#
# No refitting is involved: A* is a deterministic function of the posterior
# draws of var_between_A and the supplied V_s.

#' Sensitivity of Estimand A* to the Sampling Variances
#'
#' Recomputes Estimand A* with the supplied sampling variances multiplied by
#' each of several factors. Because A* is a deterministic function of Estimand
#' A's posterior draws and the fixed \eqn{\hat V_s}, this needs no new sampling.
#'
#' The article treats the \eqn{\hat V_s} as known and does not propagate their
#' estimation uncertainty (article Section 2.4). This function does not change
#' that assumption. It reports how sensitive the conclusion is to it. A large
#' posterior share at the zero boundary indicates a boundary-dominated,
#' nonregular A* summary; it does not make A* mathematically undefined.
#'
#' @param fit A `bhf_fit` fitted with de-attenuation enabled.
#' @param scale Numeric multipliers applied to the supplied \eqn{\hat V_s}.
#' @param prob Interval probability.
#' @return A data frame of class `bhf_astar_sensitivity`, one row per scale,
#'   with the posterior mean, interval, and boundary diagnostics for A*'s
#'   between-domain variance and proportion.
#' @seealso [variance_decomposition()], [bhf_plot_astar_sensitivity()]
#' @examples
#' \dontrun{
#' bhf_astar_sensitivity(fit)
#' bhf_astar_sensitivity(fit, scale = c(0.5, 1, 2))
#' }
#' @export
bhf_astar_sensitivity <- function(fit,
                                  scale = c(0.5, 0.75, 1, 1.25, 1.5),
                                  prob = 0.95) {
  .bhf_validate_result_fit(fit)
  if (!is.numeric(scale) || !length(scale) || anyNA(scale) ||
      any(!is.finite(scale)) || any(scale < 0)) {
    .bhf_result_abort("bhf_argument_error",
                      "`scale` must be finite, non-negative multipliers.")
  }
  interval <- .bhf_interval_contract(prob)

  sv <- fit$data$sampling_variance_info
  if (is.null(sv) || !isTRUE(sv$enabled)) {
    .bhf_result_abort(
      "bhf_result_error",
      paste0("This fit was prepared with de-attenuation disabled, so there is ",
             "no sampling-variance input to vary. Refit with ",
             "deattenuation = \"taylor\" or \"supplied\"."))
  }
  vhat <- as.numeric(fit$data$stan_data$vhat_state)
  shares <- as.numeric(fit$data$stan_data$w_state_pop_share)
  if (!length(vhat) || length(vhat) != length(shares) || anyNA(vhat) ||
      anyNA(shares) || any(!is.finite(vhat)) || any(!is.finite(shares)) ||
      any(vhat < 0) || any(shares <= 0)) {
    .bhf_result_abort(
      "bhf_fit_schema_error",
      "The fit has invalid sampling variances or population shares."
    )
  }
  base_correction <- sum(shares * vhat)

  draws <- .bhf_extract_required(
    fit, c("var_between_A", "var_within_A")
  )
  between_A <- .bhf_draw_vector(draws$var_between_A, "var_between_A")
  within_A <- .bhf_draw_vector(draws$var_within_A, "var_within_A")
  if (length(between_A) != length(within_A)) {
    .bhf_result_abort(
      "bhf_draw_shape_error",
      "A between- and within-variance fields have inconsistent draw counts."
    )
  }

  rows <- lapply(scale, function(s) {
    corr <- s * base_correction
    between <- pmax(0, between_A - corr)
    total <- between + within_A
    proportion <- ifelse(total > 0, between / total, 0)
    q <- function(v, p) unname(stats::quantile(v, p))
    lo <- unname(interval$quantiles[["lower"]])
    hi <- unname(interval$quantiles[["upper"]])
    data.frame(
      scale = s,
      correction = corr,
      between_mean = mean(between),
      between_lower = q(between, lo),
      between_median = stats::median(between),
      between_upper = q(between, hi),
      proportion_mean = mean(proportion),
      proportion_lower = q(proportion, lo),
      proportion_median = stats::median(proportion),
      proportion_upper = q(proportion, hi),
      share_at_boundary = mean(between <= 0),
      prob = interval$prob,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  attr(out, "base_correction") <- base_correction
  attr(out, "vhat_source") <- sv$provenance$mode
  class(out) <- c("bhf_astar_sensitivity", "data.frame")
  out
}

#' @export
print.bhf_astar_sensitivity <- function(x, ...) {
  cat("Estimand A* sensitivity to the supplied sampling variances\n")
  cat(sprintf("  vhat source: %s\n", attr(x, "vhat_source")))
  cat(sprintf("  baseline weighted correction: %.6g\n\n",
              attr(x, "base_correction")))
  d <- as.data.frame(x)
  show <- data.frame(
    scale = d$scale,
    `A* between` = sprintf("%.5f [%.5f, %.5f]", d$between_mean,
                           d$between_lower, d$between_upper),
    `A* proportion` = sprintf("%.4f [%.4f, %.4f]", d$proportion_mean,
                              d$proportion_lower, d$proportion_upper),
    `at zero` = sprintf("%.2f", d$share_at_boundary),
    check.names = FALSE
  )
  print(show, row.names = FALSE)
  cat(sprintf("\nIntervals are %.0f%%. 'at zero' is the posterior share of draws\n",
              100 * d$prob[1]))
  cat("truncated at the zero boundary.\n")
  invisible(x)
}

#' Plot the A* Sampling-Variance Sensitivity
#'
#' @param x The result of [bhf_astar_sensitivity()].
#' @param what Which quantity to plot.
#' @return A `ggplot` object.
#' @seealso [bhf_astar_sensitivity()]
#' @examples
#' \dontrun{
#' bhf_plot_astar_sensitivity(bhf_astar_sensitivity(fit))
#' }
#' @export
bhf_plot_astar_sensitivity <- function(x, what = c("proportion", "between")) {
  .bhf_need_ggplot2("bhf_plot_astar_sensitivity")
  what <- match.arg(what)
  if (!inherits(x, "bhf_astar_sensitivity")) {
    .bhf_result_abort("bhf_argument_error",
                      "`x` must come from bhf_astar_sensitivity().")
  }
  d <- as.data.frame(x)
  d$y  <- d[[paste0(what, "_mean")]]
  d$lo <- d[[paste0(what, "_lower")]]
  d$hi <- d[[paste0(what, "_upper")]]
  lab <- if (what == "proportion") "A* proportion between domains" else
    "A* between-domain variance"

  ggplot2::ggplot(d, ggplot2::aes(x = .data$scale, y = .data$y)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                         alpha = 0.18) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
    ggplot2::labs(
      x = "Multiplier applied to the supplied sampling variances",
      y = lab, title = "Sensitivity of A* to the fixed-variance assumption",
      caption = paste("The dashed line is the value actually used.",
                      "A* conditions on these variances as known.")) +
    ggplot2::theme_bw()
}
