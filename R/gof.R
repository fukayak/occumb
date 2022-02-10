#' @title Goodness-of-fit assessment.
#' 
#' @description \code{gof()} computes some fit statistics and their Bayesian \eqn{p}-values from the model fitting results and the original data.
#' @details
#'  TODO: Add details.
#'      - This function conducts the procedure of gof assessment described in (ref).
#'      - Two statistics: deviance and Freeman_Tukey statistics.
#'      - Bayesian p-value: I(obs < rep) / N.
#'          Extreme values may indicate inadequate model fit.
#'      - References for gof stats and Bayesian p-values.
#' @param fit An \code{occumbFit} object.
#' @param data An \code{occumbData} object used for the model fitting.
#' @param plot Logical, determine if draw scatter plots of the fit statistics.
#' @section References:
#'  TODO: Add references.
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
#' gof_result <- gof(fit, data)
#' }
#' @export
gof <- function(fit, data, plot = TRUE) {

    # Validate arguments
    qc_gof(fit, data)

    # Set constants
    y <- get_data(data, "y")
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
                dev_obs_m[j, k] <- -2 * loglik(y[, j, k], N[j, k], pi[m, , j, k])
                dev_rep_m[j, k] <- -2 * loglik(y_rep[m, , j, k], N[j, k], pi[m, , j, k])
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

# Validation for the inputs
qc_gof <- function(fit, data) {
    # Check object classes
    if (!inherits(fit, "occumbFit"))
        stop("An occumbFit class object is expected for fit")
    if (!inherits(data, "occumbData"))
        stop("An occumbData class object is expected for data")

    # Check data dimensions
    y  <- get_data(data, "y")
    pi <- get_post_samples(fit, "pi")
    if (!identical(dim(y), dim(pi)[-1]))
        stop("Dimension mismatch between the data and posterior: make sure to supply the data to which the model was applied")
}

# Fit statistics --------------------------------------------------------------
Freeman_Tukey <- function(y, N, pi) {
    sum((sqrt(y) - sqrt(N * pi))^2)
}

loglik <- function(y, N, pi) {
    stats::dmultinom(y, N, pi, log = TRUE)
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

