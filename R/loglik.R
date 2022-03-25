#' @title Pointwise log-likelihood of the fitted model.
#' 
#' @description \code{loglik()} extracts the pointwise log-likelihood matrix from
#'  an \code{occumbFit} object.
#' @details
#'  The pointwise log-likelihood is the log-likelihood of each data point
#'  evaluated with the posterior samples of parameter values. It can be used to
#'  obtain criteria for model comparisons, including WAIC and PSIS-LOO. These
#'  criteria can be conveniently calculated by applying functions from the
#'  \href{https://cran.r-project.org/web/packages/loo/index.html}{\code{loo}} package
#'  to the output of this function: see Vehtari et al. (2017) and the
#'  documentation of the \code{loo} package for details.
#' @param fit An \code{occumbFit} object.
#' @return An \eqn{M \times L} matrix, where \eqn{M} is the size of the
#'  posterior sample and \eqn{L} is the number of replicates.
#' @section References:
#'  A. Vehtari, A. Gelman, J. Gabry (2017) Practical Bayesian model evaluation
#'  using leave-one-out cross-validation and WAIC.
#'  \emph{Statistics and Computing}, \strong{27}, 1413--1432.
#'  \url{https://doi.org/10.1007/s11222-016-9696-4}
#' @examples
#' \dontrun{
#' # Generate the smallest random dataset (2 species * 2 sites * 2 reps)
#' I <- 2 # Number of species
#' J <- 2 # Number of sites
#' K <- 2 # Number of replicates
#' data <- occumbData(
#'     y = array(sample.int(I * J * K), dim = c(I, J, K)),
#'     spec_cov = list(cov1 = rnorm(I)),
#'     site_cov = list(cov2 = rnorm(J),
#'                     cov3 = factor(1:J)),
#'     repl_cov = list(cov4 = matrix(rnorm(J * K), J, K)))
#' 
#' # Fitting a null model
#' fit <- occumb(data = data)
#' 
#' # Extract log-likelihood
#' ll <- loglik(fit)  # could be passed to loo() or waic() in the loo package
#' }
#' @export
loglik <- function(fit) {

    # Validate arguments
    qc_occumbFit(fit)

    # Set constants
    y <- get_data(fit, "y")
    I <- dim(y)[1]; J <- dim(y)[2]; K <- dim(y)[3]
    N <- apply(y, c(2, 3), sum)

    pi <- get_post_samples(fit, "pi")
    M  <- dim(pi)[1]

    # Calculate fit statistics
    out <- matrix(nrow = M, ncol = J * K)
    for (m in seq_len(M)) {
        ll <- matrix(nrow = J, ncol = K)
        for (j in seq_len(J)) {
            for (k in seq_len(K)) {
                ll[j, k] <- llmulti(y[, j, k], N[j, k], pi[m, , j, k])
            }
        }
        out[m, ] <- c(ll)
    }

    # Output
    out
}

