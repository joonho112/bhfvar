# Plotting helpers. ggplot2 is a Suggests dependency: every entry point checks
# for it and fails with an actionable message rather than a missing-object error.

.bhf_ggplot2_available <- function() {
  requireNamespace("ggplot2", quietly = TRUE)
}

.bhf_need_ggplot2 <- function(fun) {
  if (!.bhf_ggplot2_available()) {
    .bhf_result_abort(
      "bhf_missing_suggest",
      sprintf(paste0(
        "%s() needs the 'ggplot2' package, which is a suggested dependency.\n",
        "  * install it with install.packages(\"ggplot2\"), or\n",
        "  * build the plot yourself from the tidy extractor output."), fun)
    )
  }
}

.bhf_as_decomposition <- function(x, prob) {
  if (inherits(x, "bhf_variance_decomposition")) return(x)
  if (inherits(x, "bhf_fit")) return(variance_decomposition(x, prob = prob))
  .bhf_result_abort(
    "bhf_argument_error",
    "`x` must be a bhf_fit or the result of variance_decomposition()."
  )
}

.bhf_interval_caption <- function(prob) {
  sprintf("Points are posterior means; bars are %.0f%% intervals.", 100 * prob)
}

#' Plot the Probability-Scale Variance Decomposition
#'
#' Shows the between-domain, within-domain, and total variance for Estimands
#' A, A*, and B side by side, so the effect of the design (A versus B) and of
#' the de-attenuation correction (A versus A*) are visible at once.
#'
#' @param x A `bhf_fit`, or the result of [variance_decomposition()].
#' @param components Which variance components to show.
#' @param prob Interval probability, used when `x` is a `bhf_fit`.
#' @return A `ggplot` object.
#' @seealso [variance_decomposition()], [bhf_plot_icc()]
#' @examples
#' \dontrun{
#' fit <- bhf_fit(prepared, model = model)
#' bhf_plot_variance(fit)
#' }
#' @export
bhf_plot_variance <- function(x,
                              components = c("between", "within", "total"),
                              prob = 0.95) {
  .bhf_need_ggplot2("bhf_plot_variance")
  components <- match.arg(components, several.ok = TRUE)
  vd <- .bhf_as_decomposition(x, prob)

  # The summary table also carries a "gaps" section whose rows reuse the
  # component names; restrict to the three estimands so those do not leak in
  # as an NA category with negative values.
  estimands <- c("A", "A_star", "B")
  tab <- vd$summary_table
  tab <- tab[tab$scale == "probability" &
               tab$estimand %in% estimands &
               tab$component %in% components &
               !is.na(tab$available) & tab$available, ]
  if (!nrow(tab)) {
    .bhf_result_abort("bhf_result_error",
                      "No probability-scale components available to plot.")
  }
  tab$estimand <- factor(tab$estimand, levels = estimands,
                         labels = c("A", "A*", "B"))
  tab$component <- factor(tab$component, levels = components)

  ggplot2::ggplot(tab, ggplot2::aes(x = .data$estimand, y = .data$mean)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper), width = 0.16) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::facet_wrap(~ .data$component, scales = "free_y") +
    ggplot2::labs(
      x = "Estimand", y = "Variance (probability scale)",
      title = "Probability-scale variance decomposition",
      caption = .bhf_interval_caption(tab$prob[1])) +
    ggplot2::theme_bw()
}

#' Plot Latent-Scale Intraclass Correlations
#'
#' A latent-scale diagnostic formed from the model's pre-centering random-effect
#' scale parameters and the logistic residual variance. These shares are not
#' empirical variances of the centered realized effects and do not by themselves
#' determine whether a random-effect term is warranted.
#'
#' @param x A `bhf_fit`, or the result of [variance_decomposition()].
#' @param prob Interval probability, used when `x` is a `bhf_fit`.
#' @return A `ggplot` object.
#' @seealso [variance_decomposition()], [bhf_plot_variance()]
#' @examples
#' \dontrun{
#' bhf_plot_icc(fit)
#' }
#' @export
bhf_plot_icc <- function(x, prob = 0.95) {
  .bhf_need_ggplot2("bhf_plot_icc")
  vd <- .bhf_as_decomposition(x, prob)

  tab <- vd$summary_table
  tab <- tab[tab$scale == "latent" & grepl("^icc_", tab$component), ]
  if (!nrow(tab)) {
    .bhf_result_abort("bhf_result_error", "No latent ICC components available.")
  }
  lvl <- c("icc_state", "icc_stratum", "icc_psu")
  tab <- tab[tab$component %in% lvl, ]
  tab$component <- factor(tab$component, levels = rev(lvl),
                          labels = rev(c("Domain", "Stratum", "PSU")))

  ggplot2::ggplot(tab, ggplot2::aes(x = .data$component, y = .data$mean)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper), width = 0.16) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL, y = "Intraclass correlation (logit scale)",
      title = "Latent scale-parameter shares",
      caption = .bhf_interval_caption(tab$prob[1])) +
    ggplot2::theme_bw()
}

#' Plot Domain Estimates
#'
#' A caterpillar plot of the domain-level probabilities with their posterior
#' intervals, ordered by the posterior mean.
#'
#' @param fit A `bhf_fit`.
#' @param estimand `"A"` or `"B"`.
#' @param prob Interval probability.
#' @param n_domains Optional cap on how many domains to draw, taking the
#'   extremes at each end. `NULL` draws all of them. For an odd cap, the extra
#'   domain is taken from the upper end of the ordered estimates.
#' @return A `ggplot` object.
#' @seealso [domain_estimates()], [bhf_plot_shrinkage()]
#' @examples
#' \dontrun{
#' bhf_plot_domains(fit, estimand = "B")
#' }
#' @export
bhf_plot_domains <- function(fit, estimand = c("A", "B"),
                             prob = 0.95, n_domains = NULL) {
  .bhf_need_ggplot2("bhf_plot_domains")
  estimand <- .bhf_match_choice(estimand, c("A", "B"), "estimand")
  .bhf_validate_result_fit(fit)
  d <- as.data.frame(domain_estimates(fit, estimand = estimand, prob = prob))
  d <- d[order(d$mean), ]

  if (!is.null(n_domains)) {
    if (!is.numeric(n_domains) || length(n_domains) != 1L ||
        is.na(n_domains) || !is.finite(n_domains) || n_domains < 2 ||
        n_domains != floor(n_domains)) {
      .bhf_result_abort("bhf_argument_error",
                        paste0("`n_domains` must be one finite integer of ",
                               "at least 2."))
    }
    k <- min(as.integer(n_domains), nrow(d))
    n_low <- k %/% 2L
    n_high <- k - n_low
    keep <- c(seq_len(n_low),
              seq.int(nrow(d) - n_high + 1L, nrow(d)))
    d <- d[keep, ]
  }
  d$domain <- factor(d$domain, levels = d$domain)

  ggplot2::ggplot(d, ggplot2::aes(x = .data$domain, y = .data$mean)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper), width = 0) +
    ggplot2::geom_point(size = 1.9) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL, y = sprintf("Estimand %s (probability)", estimand),
      title = "Domain estimates",
      caption = .bhf_interval_caption(prob)) +
    ggplot2::theme_bw()
}

#' Plot Shrinkage of Domain Estimates
#'
#' Compares the weighted raw proportion in each domain against the model-based
#' estimate and its posterior interval. Differences may reflect partial pooling
#' as well as standardization over the fitted survey-design effects.
#'
#' @param fit A `bhf_fit`.
#' @param estimand `"A"` or `"B"`.
#' @param prob Interval probability.
#' @return A `ggplot` object.
#' @seealso [domain_estimates()], [bhf_plot_domains()]
#' @examples
#' \dontrun{
#' bhf_plot_shrinkage(fit)
#' }
#' @export
bhf_plot_shrinkage <- function(fit, estimand = c("A", "B"),
                               prob = 0.95) {
  .bhf_need_ggplot2("bhf_plot_shrinkage")
  estimand <- .bhf_match_choice(estimand, c("A", "B"), "estimand")
  .bhf_validate_result_fit(fit)
  ad <- fit$data$analysis_data
  if (is.null(ad) || !all(c("y", "raw_weight", "state_id") %in% names(ad))) {
    .bhf_result_abort(
      "bhf_result_error",
      "The fit does not carry the analysis rows needed for a shrinkage plot.")
  }

  # domain_estimates() already carries `n`; do not duplicate the column here or
  # the merge below renames both to n.x / n.y.
  raw <- do.call(rbind, lapply(sort(unique(ad$state_id)), function(s) {
    z <- ad[ad$state_id == s, ]
    data.frame(domain_id = s,
               raw = sum(z$raw_weight * z$y) / sum(z$raw_weight),
               stringsAsFactors = FALSE)
  }))
  d <- merge(as.data.frame(domain_estimates(fit, estimand = estimand, prob = prob)),
             raw, by = "domain_id")

  ggplot2::ggplot(d, ggplot2::aes(x = .data$raw, y = .data$mean)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey55") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper), width = 0) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n), alpha = 0.75) +
    ggplot2::scale_size_continuous(name = "Domain n") +
    ggplot2::labs(
      x = "Weighted raw proportion", y = "Model-based estimate",
      title = "Model-based versus weighted direct domain estimates",
      caption = paste(
        .bhf_interval_caption(prob),
        "Differences can reflect partial pooling and design-effect standardization."
      )) +
    ggplot2::theme_bw()
}
