#' bhfvar: Bayesian Hybrid Framework for Variance Decomposition
#'
#' The bhfvar package implements the Bayesian Hybrid Framework for variance
#' decomposition in complex surveys with post-hoc domains. It provides tools
#' for separating substantive geographic variation from design artifacts and
#' sampling noise.
#'
#' @section Key Features:
#' \itemize{
#'   \item Bayesian Pseudo-Likelihood estimation for design consistency
#'   \item Hybrid generalized linear mixed models with domain and PSU effects
#'   \item Dual Estimand Framework: Policy (A/A*) and Descriptive (B) estimands
#'   \item De-attenuation for finite-sample variance inflation correction
#'   \item Comprehensive diagnostic and visualization tools
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
#' @section Design Philosophy:
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
