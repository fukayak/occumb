#' @title Goodness-of-fit assessment of the fitted model.
#' 
#' @description \code{gof()} calculates some omnibus discrepancy measures and 
#'  their Bayesian \eqn{p}-values for the fitted model using the posterior
#'  predictive check approach.
#' @details
#'  A discrepancy statistic of the fitted model is obtained using the procedure
#'  of posterior predictive checking.
#'  The following statistics are currently available:
#'      \describe{
#'          \item{Freeman-Tukey statistics (default)}{\eqn{T(\boldsymbol{y}, \boldsymbol{\theta}) = \sum_{i}\sum_{j}\sum_{k}\left(\sqrt{y_{ijk}} - \sqrt{E(y_{ijk} \mid \boldsymbol{\theta})}\right)^2}{\sum_i \sum_j \sum_k (sqrt(y[i, j, k]) - sqrt(E(y[i, j, k] | theta)))^2}}
#'          \item{Deviance statistics}{\eqn{T(\boldsymbol{y}, \boldsymbol{\theta}) = -2 \log p(\boldsymbol{y} \mid \boldsymbol{\theta})}{T(y, theta) = -2 * log(p(y | theta))}}
#'      }
#'  where \eqn{\boldsymbol{y} = \{y_{ijk}\}}{y}, \eqn{\boldsymbol{\theta}}{theta}, 
#'  \eqn{E(y_{ijk} \mid \boldsymbol{\theta})}{E(y[i, j, k] | theta)}, and
#'  \eqn{\log p(\boldsymbol{y} \mid \boldsymbol{\theta})}{log(p(y | theta))}
#'  are sequence read count data, parameters and latent variables of the model,
#'  expected value of the data conditional on \eqn{\boldsymbol{\theta}}{theta},
#'  and log-likelihood of the model, respectively.
#'  The Bayesian \eqn{p}-value is estimated as the probability that the value of
#'  the discrepancy statistics of replicated dataset is more extreme than that
#'  of the observed dataset.
#'  An extreme Bayesian \eqn{p}-value may indicate an inadequate model fit.
#'  See, e.g., Gelman et al. (2014), Kéry and Royle (2016), and
#'  Conn et al. (2018) for more details on the procedures for posterior
#'  predictive checking.
#' @param fit An \code{occumbFit} object.
#' @param stats The discrepancy statistics to be applied.
#' @param plot Logical, determine if draw scatter plots of the fit statistics.
#' @return A list with the following named elements:
#'      \describe{
#'          \item{\code{stats}}{The discrepancy statistics applied.}
#'          \item{\code{p_value}}{Bayesian \eqn{p}-value.}
#'          \item{\code{stats_obs}}{Discrepancy statistics for observed data.}
#'          \item{\code{stats_rep}}{Discrepancy statistics for repeated data.}
#'      }
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
gof <- function(fit,
                stats = c("Freeman_Tukey", "deviance"),
                plot = TRUE) {

    # Validate arguments
    qc_occumbFit(fit)
    stats <- match.arg(stats)

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
    stats_obs <- stats_rep <- vector(length = M)
    for (m in seq_len(M)) {
        stats_obs_m <- stats_rep_m <- matrix(nrow = J, ncol = K)
        for (j in seq_len(J)) {
            for (k in seq_len(K)) {
                if (stats == "Freeman_Tukey") {
                    stats_obs_m[j, k]  <- Freeman_Tukey(y[, j, k], N[j, k], pi[m, , j, k])
                    stats_rep_m[j, k]  <- Freeman_Tukey(y_rep[m, , j, k], N[j, k], pi[m, , j, k])
                }
                if (stats == "deviance") {
                    stats_obs_m[j, k] <- -2 * llmulti(y[, j, k], N[j, k], pi[m, , j, k])
                    stats_rep_m[j, k] <- -2 * llmulti(y_rep[m, , j, k], N[j, k], pi[m, , j, k])
                }
            }
        }
        stats_obs[m] <- sum(stats_obs_m)
        stats_rep[m] <- sum(stats_rep_m)
    }

    # Output (plot and object)
    if (plot) plot_gof(stats_obs, stats_rep, stats)
    out <- list(stats = stats,
                p_value = Bayesian_p_value(stats_obs, stats_rep),
                stats_obs = stats_obs,
                stats_rep = stats_rep)
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

