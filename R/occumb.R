#' @include classes.R
NULL

# Class for model-fit results of occumb
setClass("occumbFit", slots = c(fit = "jagsUI"))

#' Model-fitting function.
#' 
#' \code{occumb()} fits a specified site-occupancy model to a dataset and returns
#' a model-fit object containing posterior samples.
#' 
#' @param phi_formula A right-hand side formula describing covariates of
#'        sequence relative dominance.
#' @param theta_formula A right-hand side formula describing covariates of
#'        sequence capture probability.
#' @param psi_formula A right-hand side formula describing covariates of
#'        occupancy probability.
#' @param data occumbData class object.
#' @param prior_prec Precision of the normal prior distribution for the
#'        community-level average of species-specific parameters.
#' @param prior_ulim Upper limit of the uniform prior distribution for the
#'        standard deviation of species-specific parameters.
#' @param n.chains Number of Markov chains to run.
#'        Passed to the \code{n.chains} argument of jagsUI's \code{jags()}
#'        function.
#' @param n.adapt Number of iterations to run in the JAGS adaptive phase.
#'        Passed to the \code{n.adapt} argument of jagsUI's \code{jags()}
#'        function.
#'        See the document of jagsUI's \code{jags()} function for details.
#' @param n.burnin Number of iterations at the beginning of the chain to discard.
#'        Passed to the \code{n.burnin} argument of jagsUI's \code{jags()}
#'        function.
#'        See the document of jagsUI's \code{jags()} function for details.
#' @param n.thin Thinning rate. Must be a positive integer.
#'        Passed to the \code{n.thin} argument of jagsUI's \code{jags()}
#'        function.
#' @param n.iter Total number of iterations per chain (including burn-in).
#'        Passed to the \code{n.iter} argument of jagsUI's \code{jags()}
#'        function.
#' @param parallel If TRUE, run MCMC chains in parallel on multiple CPU cores.
#'        Passed to the \code{parallel} argument of jagsUI's \code{jags()}
#'        function.
#' @section Details:
#'      \code{occumb()} allows the fitting of a range of site occupancy models
#'      for sequence read counts that includes covariates at different levels
#'      of the data generation process, via the specification of the model
#'      formula in \code{phi_formula}, \code{theta_formula}, and
#'      \code{psi_formula.}
#'      The specified model and dataset are passed to JAGS via the
#'      interface provided by the jagsUI package, where Markov chain
#'      Monte Carlo methods are used to obtain posterior samples of model
#'      parameters.
#'
#'      For the resulting object, the following methods and functions provided
#'      by the jagsUI package can be applied:
#'
#'      plot(); print(); summary()
#'
#'      Currently, occumb() supports covariate modeling only for
#'      psi_formula in which site covariates can be incorporated.
#' @return  An S4 object of the \code{occumbFit} class containing the results of
#'          model fitting.
#' @export
occumb <- function(phi_formula = ~ 1,
                   theta_formula = ~ 1,
                   psi_formula = ~ 1,
                   data,
                   prior_prec = 1E-3,
                   prior_ulim = 1E3,
                   n.chains = 6,
                   n.adapt = 1000,
                   n.burnin = 30000,
                   n.thin = 500,
                   n.iter = 500 * n.thin + n.burnin,
                   parallel = TRUE) {
    # QC

    # Set constants
    const <- set_const(data)

    # Define arguments for the specified model
    margs <- set_modargs(phi_formula, theta_formula, psi_formula, data)

    # Set data list
    dat <- list(I          = const$I,
                J          = const$J,
                K          = const$K,
                N          = const$N,
                y          = const$y,
                cov_phi    = margs$cov_phi,
                cov_theta  = margs$cov_theta,
                cov_psi    = margs$cov_psi,
                M          = margs$M,
                m_phi      = margs$m_phi,
                m_theta    = margs$m_theta,
                m_psi      = margs$m_psi,
                prior_prec = prior_prec,
                prior_ulim = prior_ulim)

    # Set initial values
    inits <- function() {
        list(z = matrix(1, const$I, const$J),
             u = array(1, dim = c(const$I, const$J, const$K)),
             x = array(stats::rnorm(const$I * const$J * const$K,
                                    mean = 1, sd = 0.1),
                       dim = c(const$I, const$J, const$K)),
             spec_eff = matrix(stats::rnorm(const$I * margs$M, sd = 0.1),
                               const$I, margs$M),
             Mu       = stats::rnorm(margs$M, sd = 0.1),
             sigma    = stats::rnorm(margs$M, mean = 1, sd = 0.1),
             rho      = matrix(stats::rnorm(margs$M^2, sd = 0.1),
                               margs$M, margs$M))
    }

    # Set parameters monitored
    params <- c("Mu", "sigma", "rho", "alpha", "beta", "gamma",
                "phi", "theta", "psi", "z", "pi")

    # Run MCMC in JAGS
    fit <- jagsUI::jags(
        dat, inits, params,
        system.file("jags",
                    sprintf("occumb%s.jags", margs$code),
                    package = "occumb"),
        n.chains = n.chains,
        n.adapt  = n.adapt,
        n.iter   = n.iter,
        n.burnin = n.burnin,
        n.thin   = n.thin,
        parallel = parallel)

    # Output
    out <- methods::new("occumbFit", fit = fit)
    out
}

# Extract constants for the model
set_const <- function(data) {
    y <- data@y
    I <- dim(y)[1]              # Number of species
    J <- dim(y)[2]              # Number of sites
    K <- dim(y)[3]              # Number of replicates
    N <- apply(y, c(2, 3), sum) # Sequence depth

    for (j in 1:J) {
        for (k in 1:K) {
            # Set y to NA and N to 1 when sequence depth is zero
            if (N[j, k] == 0) {
                y[, j, k] <- NA
                N[j, k]   <- 1
            }
        }
    }

    out <- list(y = y, I = I, J = J, K = K, N = N)
    out
}

# Set model-specific arguments
set_modargs <- function(phi_formula, theta_formula, psi_formula, data) {
    M    <- 0           # Order of species effects
    code <- "000000"    # Model code

    if (phi_formula == ~ 1) {
        cov_phi <- 1
        M       <- M + 1
        m_phi   <- 1
        substr(code, 1, 1) <- "1"
    } else {
        stop("Covariates in phi_formula are not yet supported, sorry.")
    }

    if (theta_formula == ~ 1) {
        cov_theta <- 1
        M         <- M + 1
        m_theta   <- seq(m_phi + 1, m_phi + 1)
        substr(code, 2, 2) <- "1"
    } else {
        stop("Covariates in theta_formula are not yet supported, sorry.")
    }

    if (psi_formula == ~ 1) {
        cov_psi <- 1
        M       <- M + 1
        m_psi   <- seq(m_theta + 1, m_theta + 1)
        substr(code, 3, 3) <- "1"
    } else {
        cov_psi <- stats::model.matrix(stats::as.formula(psi_formula),
                                       data@site_cov)
        M       <- M + ncol(cov_psi)
        m_psi   <- seq(m_theta + 1, m_theta + ncol(cov_psi))
        substr(code, 3, 3) <- "2"
    }

    out <- list(M         = M,
                code      = code,
                cov_phi   = cov_phi,
                cov_theta = cov_theta,
                cov_psi   = cov_psi,
                m_phi     = m_phi,
                m_theta   = m_theta,
                m_psi     = m_psi)
    out
}

# Auto-generate JAGS model code
write_jags_model <- function(phi, theta, psi, shared_effects) {

    model <- readLines(system.file("jags",
                                   "occumb_template1.jags",
                                   package = "occumb"))

    if (phi == "i")
        model <- c(model,
                   "                x[i, j, k] ~ dgamma(phi[i], 1)")
    if (phi == "ij")
        model <- c(model,
                   "                x[i, j, k] ~ dgamma(phi[i, j], 1)")
    if (phi == "ijk")
        model <- c(model,
                   "                x[i, j, k] ~ dgamma(phi[i, j, k], 1)")

    model <- c(model,
               readLines(system.file("jags",
                                     "occumb_template2.jags",
                                     package = "occumb")))

    if (theta == "i")
        model <- c(model,
                   "                u[i, j, k] ~ dbern(z[i, j] * theta[i])")
    if (theta == "ij")
        model <- c(model,
                   "                u[i, j, k] ~ dbern(z[i, j] * theta[i, j])")
    if (theta == "ijk")
        model <- c(model,
                   "                u[i, j, k] ~ dbern(z[i, j] * theta[i, j, k])")

    model <- c(model,
               readLines(system.file("jags",
                                     "occumb_template3.jags",
                                     package = "occumb")))

    if (psi == "i")
        model <- c(model,
                   "            z[i, j] ~ dbern(psi[i])")
    if (psi == "ij")
        model <- c(model,
                   "            z[i, j] ~ dbern(psi[i, j])")

    model <- c(model,
               readLines(system.file("jags",
                                     "occumb_template4.jags",
                                     package = "occumb")))

    if (phi == "i")
        model <- c(model,
                   "        log(phi[i])     <- inprod(alpha[i, ], cov_phi[])")
    if (phi == "ij")
        model <- c(model,
                   "        log(phi[i, j])     <- inprod(alpha[i, ], cov_phi[i, j, ])")
    if (phi == "ijk")
        model <- c(model,
                   "        log(phi[i, j, k])     <- inprod(alpha[i, ], cov_phi[i, j, k, ])")

    if (theta == "i")
        model <- c(model,
                   "        logit(theta[i]) <- inprod(beta[i, ], cov_theta[])")
    if (theta == "ij")
        model <- c(model,
                   "        logit(theta[i, j]) <- inprod(beta[i, ], cov_theta[i, j, ])")
    if (theta == "ijk")
        model <- c(model,
                   "        logit(theta[i, j, k]) <- inprod(beta[i, ], cov_theta[i, j, k, ])")

    if (psi == "i")
        model <- c(model,
                   "        logit(psi[i])   <- inprod(gamma[i, ], cov_psi[])")
    if (psi == "ij")
        model <- c(model,
                   "        logit(psi[i, j])   <- inprod(gamma[i, ], cov_psi[i, j, ])")

    model <- c(model,
               readLines(system.file("jags",
                                     sprintf("occumb_template5%s.jags",
                                             ifelse(shared_effects, "b", "a")),
                                     package = "occumb")))

    model
}

