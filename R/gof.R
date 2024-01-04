# Data format class for occumb
setClass("occumbGof",
         slots = c(stats = "character",
                   p_value = "numeric",
                   stats_obs = "numeric",
                   stats_rep = "numeric"))

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
#'          \item{Freeman-Tukey statistics (default)}{\eqn{T_{\textrm{FT}} = \sum_{i,j,k}\left(\sqrt{y_{ijk}} - \sqrt{E(y_{ijk} \mid \pi_{ijk})}\right)^2}{T_{FT} = \sum_{i, j, k} (sqrt(y[i, j, k]) - sqrt(E(y[i, j, k] | pi[i, j, k])))^2}}
#'          \item{Deviance statistics}{\eqn{T_{\textrm{deviance}} = -2 \sum_{j,k} \log \textrm{Multinomial}(\boldsymbol{y}_{jk} \mid \boldsymbol{\pi}_{jk})}{T_{deviance} = -2 * \sum_{j, k} log(Multinomial(y[, j, k] | pi[, j, k]))}}
#'      }
#'  where \eqn{i}, \eqn{j}, and \eqn{k} are the subscripts of species, site, and
#'  replicate, respectively,
#'  \eqn{y_{ijk}}{y} is sequence read count data,
#'  \eqn{\pi_{ijk}}{pi} is multinomial cell probabilities of sequence read
#'  counts,
#'  \eqn{E(y_{ijk} \mid \pi_{ijk})}{E(y[i, j, k] | pi[i, j, k])} is expected
#'  value of the sequence read counts conditional on their cell probabilities,
#'  and \eqn{\log \textrm{Multinomial}(\boldsymbol{y}_{jk} \mid \boldsymbol{\pi}_{jk})}{log(Multinomial(y[, j, k] | pi[, j, k]))}
#'  is the multinomial log-likelihood of the sequence read counts in replicate
#'  \eqn{k}{k} of site \eqn{j}{j} conditional on their cell probabilities.
#'
#'  The Bayesian \eqn{p}-value is estimated as the probability that the value of
#'  the discrepancy statistics of replicated dataset is more extreme than that
#'  of the observed dataset.
#'  An extreme Bayesian \eqn{p}-value may indicate an inadequate model fit.
#'  See, e.g., Gelman et al. (2014), Kéry and Royle (2016), and
#'  Conn et al. (2018) for more details on the procedures for posterior
#'  predictive checking.
#'
#'  Computations can be run in parallel on multiple CPU cores where the `cores`
#'  argument controls the degree of parallelization.
#' @param fit An \code{occumbFit} object.
#' @param stats The discrepancy statistics to be applied.
#' @param cores The number of cores to use for parallelization.
#' @param plot Logical, determine if draw scatter plots of the fit statistics.
#' @param ... Additional arguments passed to the default \code{plot} method.
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
#'  \doi{10.1002/ecm.1314}
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
#' \donttest{
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
#' gof_result
#' }
#' @export
gof <- function(fit,
                stats = c("Freeman_Tukey", "deviance"),
                cores = 1L,
                plot = TRUE, ...) {

    # Validate arguments
    assert_occumbFit(fit)
    stats <- match.arg(stats)

    # Set constants
    y <- get_data(fit, "y")
    I <- dim(y)[1]; J <- dim(y)[2]; K <- dim(y)[3]
    N <- apply(y, c(2, 3), sum)

    pi <- get_post_samples(fit, "pi")
    M  <- dim(pi)[1]

    # Generate replicate data and calculate fit statistics
    if (cores == 1) {
        y_rep <- lapply(X = seq_len(M),
                        FUN = get_y_rep,
                        y = y,
                        N = N,
                        pi = pi)
        stats_obs <- unlist(
            lapply(X = seq_len(M),
                   FUN = get_stats,
                   stats = stats,
                   y = y,
                   N = N,
                   pi = pi)
        )
        stats_rep <- unlist(
            lapply(X = seq_len(M),
                   FUN = .get_stats,
                   stats = stats,
                   y_rep = y_rep,
                   N = N,
                   pi = pi)
        )
    } else {
        if (.Platform$OS.type == "windows") {
            # On Windows use makePSOCKcluster() and parLapply() for multiple cores
            cl <- parallel::makePSOCKcluster(cores)
            parallel::clusterEvalQ(cl, library(occumb))
            on.exit(parallel::stopCluster(cl))
            y_rep <- parallel::parLapply(cl = cl,
                                         X = seq_len(M),
                                         fun = get_y_rep,
                                         y = y,
                                         N = N,
                                         pi = pi)

            stats_obs <- unlist(
                parallel::parLapply(cl = cl,
                                    X = seq_len(M),
                                    fun = get_stats,
                                    stats = stats,
                                    y = y,
                                    N = N,
                                    pi = pi)
            )

            stats_rep <- unlist(
                parallel::parLapply(cl = cl,
                                    X = seq_len(M),
                                    fun = .get_stats,
                                    stats = stats,
                                    y_rep = y_rep,
                                    N = N,
                                    pi = pi)
            )
        } else {
            # On Mac or Linux use mclapply() for multiple cores
            y_rep <- parallel::mclapply(mc.cores = cores,
                                        X = seq_len(M),
                                        FUN = get_y_rep,
                                        y = y,
                                        N = N,
                                        pi = pi)

            stats_obs <- unlist(
                parallel::mclapply(mc.cores = cores,
                                   X = seq_len(M),
                                   FUN = get_stats,
                                   stats = stats,
                                   y = y,
                                   N = N,
                                   pi = pi)
            )

            stats_rep <- unlist(
                parallel::mclapply(mc.cores = cores,
                                   X = seq_len(M),
                                   FUN = .get_stats,
                                   stats = stats,
                                   y_rep = y_rep,
                                   N = N,
                                   pi = pi)
            )
        }
    }

    # Output (plot and object)
    if (plot) plot_gof(stats_obs, stats_rep, stats, ...)
    out <- methods::new("occumbGof",
                        stats = stats,
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

# Generate replicate data
get_y_rep <- function(m, y, N, pi) {
    I <- dim(y)[1]; J <- dim(y)[2]; K <- dim(y)[3]
    y_rep_m <- array(dim = c(I, J, K))
    for (j in seq_len(J)) {
        for (k in seq_len(K)) {
            if (N[j, k] > 0) {
                y_rep_m[, j, k] <- stats::rmultinom(1, N[j, k], pi[m, , j, k])
            } else {
                y_rep_m[, j, k] <- 0
            }
        }
    }
    y_rep_m
}

# Calculate fit statistics
.get_stats <- function(m, stats, y_rep, N, pi) {
    get_stats(m, stats, y_rep[[m]], N, pi)
}

get_stats <- function(m, stats, y, N, pi) {
    J <- dim(y)[2]; K <- dim(y)[3]
    stats_m <- matrix(nrow = J, ncol = K)
    for (j in seq_len(J)) {
        for (k in seq_len(K)) {
            if (stats == "Freeman_Tukey") {
                stats_m[j, k]  <- Freeman_Tukey(y[, j, k], N[j, k], pi[m, , j, k])
            }
            if (stats == "deviance") {
                stats_m[j, k] <- -2 * llmulti(y[, j, k], N[j, k], pi[m, , j, k])
            }
        }
    }
    sum(stats_m)
}

# Calculate Bayesian p-values
Bayesian_p_value <- function(stat_obs, stat_rep) mean(stat_obs < stat_rep)

# Plot function
plot_gof <- function(stat_obs, stat_rep, statistics, ...) {
    pval <- Bayesian_p_value(stat_obs, stat_rep)
    stats_print <- ifelse(statistics == "Freeman_Tukey", "Freeman-Tukey", statistics)
    plot(stat_rep ~ stat_obs,
         xlim = range(c(stat_obs, stat_rep)),
         ylim = range(c(stat_obs, stat_rep)),
         main = paste(stats_print, "statistics | Bayesian p-value =", round(pval, 5)),
         xlab = paste("Statistics for observed data"),
         ylab = paste("Statistics for replicated data"), ...)
    graphics::abline(a = 0, b = 1)
}

