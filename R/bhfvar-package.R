#' bhfvar: Bayesian Hybrid Framework for Variance Decomposition
#'
#' The bhfvar package implements an article-informed crossed state/stratum
#' model with PSUs nested within strata, pseudo-likelihood survey weights, and
#' probability-scale variance decompositions for post-hoc domains.
#'
#' @section Supported workflow:
#' \itemize{
#'   \item Crossed state and stratum effects with PSU effects nested in stratum
#'   \item Globally mean-one likelihood-weight scaling (legacy D2 is deprecated)
#'   \item Probability-scale Estimands A, A*, and B with explicit gaps
#'   \item Fixed supplied or Taylor-linearized sampling variances for A*
#'   \item Versioned prepared-data, fit, and extractor result contracts
#' }
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{compile_bhf_model}}: Compile Stan model (once per session)
#'   \item \code{\link{prepare_bhf_data}}: Prepare data for Stan
#'   \item \code{\link{bhf_fit}}: Fit the BHF model
#'   \item \code{\link{variance_decomposition}}: Extract variance components
#'   \item \code{\link{domain_estimates}}: Extract domain-specific estimates
#' }
#'
#' @section Workflow:
#' The recommended workflow is:
#' \enumerate{
#'   \item Compile the Stan model once per R session using \code{compile_bhf_model()}
#'   \item Prepare your data using \code{prepare_bhf_data()}
#'   \item Fit the model using \code{bhf_fit()}
#'   \item Extract results using \code{variance_decomposition()} and \code{domain_estimates()}
#' }
#'
#' @section What has and has not been established:
#' The implementation is verified against the model specified in Appendix A of
#' the accompanying article: a frozen reference oracle, algebraic property
#' tests, and centering-invariant tests all pass within their tested contracts.
#'
#' Interval calibration is a separate question and has not been established.
#' An internal simulation study found coverage of central 90\% intervals below
#' nominal for some quantities, which is consistent with the general behaviour
#' of unadjusted pseudo-posteriors (see the interpretation boundary below).
#' Users should not assume nominal frequentist coverage. The restricted-data
#' application reported in the article is not reproduced by this package.
#'
#' @section Interpretation boundary:
#' A* conditions on supplied or Taylor-estimated sampling variances; their
#' estimation uncertainty is not propagated. Pseudo-posterior intervals and
#' pseudo log likelihood do not automatically have ordinary posterior coverage
#' or ordinary observation-level LOO interpretations.
#'
#' @section Compilation design:
#' This package uses a "defensive" programming approach where the Stan model
#' is compiled explicitly by the user once per session, rather than being
#' pre-compiled during package installation. This approach:
#' \itemize{
#'   \item Avoids rstantools caching issues
#'   \item Provides clearer error messages when compilation fails
#'   \item Ensures compatibility across different R/Stan versions
#'   \item Gives users more control over the compilation process
#' }
#'
#' @references
#' Lee, J., & Hooper, A. (2026). Disentangling signal from noise: A Bayesian
#' hybrid framework for variance decomposition in complex surveys with post-hoc
#' domains. \emph{Mathematics}, 14(3), 512. \doi{10.3390/math14030512}
#'
#' @author JoonHo Lee \email{jlee296@@ua.edu}
#'
#' @docType package
#' @name bhfvar-package
#' @aliases bhfvar
"_PACKAGE"

## usethis namespace: start
#' @importFrom rstan stan_model sampling extract
#' @importFrom posterior summarise_draws as_draws_df
#' @importFrom stats var sd quantile median
#' @importFrom methods is
## usethis namespace: end
NULL
