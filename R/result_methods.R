.bhf_print_interval <- function(interval) {
  if (!inherits(interval, "bhf_interval_contract")) return("credible interval")
  paste0(interval$interval_label, " [", interval$labels[["lower"]], ", ",
         interval$labels[["upper"]], "]")
}

#' @export
print.bhf_variance_decomposition <- function(x, ...) {
  cat("=== BHF Variance Decomposition (schema ", x$schema_version,
      ") ===\n", sep = "")
  cat("Interval: ", .bhf_print_interval(x$interval), "\n", sep = "")
  cat("Scales: latent diagnostics are separate from probability A/A*/B.\n")
  if (!isTRUE(x$A_star$available)) {
    cat("A*: unavailable -- ", x$A_star$unavailable_reason, "\n", sep = "")
  }
  display <- x$summary_table[c(
    "section", "estimand", "component", "mean", "lower", "upper",
    "available"
  )]
  print.data.frame(display, row.names = FALSE, ...)
  invisible(x)
}

#' @export
print.bhf_domain_estimates <- function(x, ...) {
  interval <- attr(x, "interval", exact = TRUE)
  cat("=== BHF Domain Estimates: Estimand ",
      attr(x, "estimand", exact = TRUE), " ===\n", sep = "")
  cat("Interval: ", .bhf_print_interval(interval), "\n", sep = "")
  display <- x
  class(display) <- "data.frame"
  print.data.frame(display, row.names = FALSE, ...)
  invisible(x)
}

#' @export
print.bhf_overall_estimate <- function(x, ...) {
  cat("=== BHF Overall Estimate: Estimand ", x$estimand, " ===\n", sep = "")
  cat("Mean: ", format(x$mean, digits = 6), "; ",
      .bhf_print_interval(x$interval), ": [",
      format(x$lower, digits = 6), ", ", format(x$upper, digits = 6),
      "]\n", sep = "")
  invisible(x)
}

#' @export
print.bhf_log_lik <- function(x, ...) {
  scope <- attr(x, "scope", exact = TRUE)
  cat("=== BHF Log-Likelihood Contributions ===\n")
  cat("Kind: ", attr(x, "kind", exact = TRUE), "; aggregation: ",
      attr(x, "aggregation", exact = TRUE), "; dimensions: ", nrow(x),
      " draws x ", ncol(x), " units\n", sep = "")
  cat("Scope: ", scope$caveat, "\n", sep = "")
  invisible(x)
}
