#' @title Goodness-of-fit assessment of the fitted model.
#' 
#' @description \code{gof()} calculates some omnibus discrepancy measures and 
#'  their Bayesian \eqn{p}-values for the fitted model using the posterior
#'  predictive check approach.
#' @details
#'  Two discrepancy measures, deviance and the Freeman-Tukey statistics, are
#'  obtained using the procedure of posterior predictive checking.
#'  The Bayesian \eqn{p}-value is estimated as the probability that the value of
#'  the discrepancy measure of replicated data is more extreme than that of the
#'  observed data.
#'  An extreme Bayesian \eqn{p}-value may indicate an inadequate model fit.
#'  See, e.g., Gelman et al. (2014), Kéry and Royle (2016), and
#'  Conn et al. (2018) for more details on the procedures for posterior
#'  predictive checking.
#' @param fit An \code{occumbFit} object.
#' @return A list with the following named elements in which results for
#'  deviance and the Freeman-Tukey statistics are recorded:
#'      \describe{
#'          \item{\code{p_values}}{Bayesian \eqn{p}-value.}
#'          \item{\code{stats_obs}}{Discrepancy measures for observed data.}
#'          \item{\code{stats_rep}}{Discrepancy measures for repeated data.}
#'      }
#' @param plot Logical, determine if draw scatter plots of the fit statistics.
#' @section References:
#'  P. B. Conn, D. S. Johnson, P. J. Williams, S. R. Melin and M. B. Hooten.
#'  (2018) A guide to Bayesian model checking for ecologists.
#'  \emph{Ecological Monographs} \strong{88}:526--542.
#'  \url{https://doi.org/10.1002/ecm.1314}
#'
#'  A. Gelman, J. B. Carlin, H. S. Stern D. B. Dunson, A. Vehtari and
#'  D. B. Rubin (2013) \emph{Bayesian Data Analysis}. 3rd edition.
#'  Chapman and Hall/CRC. \url{http://www.stat.columbia.edu/~gelman/book/}
#'
#'  M. Kéry and J. A. Royle (2016) \emph{Applied Hierarchical Modeling in
#'  Ecology --- Analysis of Distribution, Abundance and Species Richness in R
#'  and BUGS. Volume 1: Prelude and Static Models}. Academic Press.
#'  \url{https://www.mbr-pwrc.usgs.gov/pubanalysis/keryroylebook/}
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
#' # Goodness-of-fit assessment
#' gof_result <- gof(fit)
#' gof_result$p_values  # print p-values
#' }
#' @export
gof <- function(fit, plot = TRUE) {

    # Validate arguments
    qc_occumbFit(fit)

    # Set constants
    y <- get_data(fit, "y")
    I <- dim(y)[1]; J <- dim(y)[2]; K <- dim(y)[3]
    N <- apply(y, c(2, 3), sum)

    pi <- get_post_samples(fit, "pi")
    M  <- dim(pi)[1]

    # Generate replicate data
    y_rep <- array(dim = c(M, I, J, K))
    for (m in seq_len(M)) {
        for (j in seq_len(J)) {
            for (k in seq_len(K)) {
                if (N[j, k] > 0) {
                    y_rep[m, , j, k] <- stats::rmultinom(1, N[j, k], pi[m, , j, k])
                } else {
                    y_rep[m, , j, k] <- 0
                }
            }
        }
    }

    # Calculate fit statistics
    dev_obs <- dev_rep <- FT_obs <- FT_rep <- vector(length = M)
    for (m in seq_len(M)) {
        dev_obs_m <- dev_rep_m <- FT_obs_m <- FT_rep_m <- matrix(nrow = J, ncol = K)
        for (j in seq_len(J)) {
            for (k in seq_len(K)) {
                dev_obs_m[j, k] <- -2 * llmulti(y[, j, k], N[j, k], pi[m, , j, k])
                dev_rep_m[j, k] <- -2 * llmulti(y_rep[m, , j, k], N[j, k], pi[m, , j, k])
                FT_obs_m[j, k]  <- Freeman_Tukey(y[, j, k], N[j, k], pi[m, , j, k])
                FT_rep_m[j, k]  <- Freeman_Tukey(y_rep[m, , j, k], N[j, k], pi[m, , j, k])
            }
        }
        dev_obs[m] <- sum(dev_obs_m)
        dev_rep[m] <- sum(dev_rep_m)
        FT_obs[m]  <- sum(FT_obs_m)
        FT_rep[m]  <- sum(FT_rep_m)
    }

    # Output (plot and object)
    if (plot) {
        graphics::par(mfrow = c(1, 2))
        plot_gof(dev_obs, dev_rep, "Deviance")
        plot_gof(FT_obs, FT_rep, "Freeman-Tukey")
    }

    out <- list(p_values = list(
                    deviance = Bayesian_p_value(dev_obs, dev_rep),
                    Freeman_Tukey = Bayesian_p_value(FT_obs, FT_rep)),
                stats_obs = list(
                    deviance = dev_obs,
                    Freeman_Tukey = FT_obs),
                stats_rep = list(
                    deviance = dev_rep,
                    Freeman_Tukey = FT_rep))
    out
}

# Fit statistics --------------------------------------------------------------
Freeman_Tukey <- function(y, N, pi) {
    sum((sqrt(y) - sqrt(N * pi))^2)
}
# -----------------------------------------------------------------------------

# Calculate Bayesian p-values
Bayesian_p_value <- function(stat_obs, stat_rep) mean(stat_obs < stat_rep)

# Plot function
plot_gof <- function(stat_obs, stat_rep, statistics) {
    pval <- Bayesian_p_value(stat_obs, stat_rep)
    plot(stat_rep ~ stat_obs,
         xlim = range(c(stat_obs, stat_rep)),
         ylim = range(c(stat_obs, stat_rep)),
         main = paste(statistics, "| Bayesian p-value =", pval),
         xlab = paste0(statistics, "_obs"),
         ylab = paste0(statistics, "_rep"))
    graphics::abline(a = 0, b = 1)
}

