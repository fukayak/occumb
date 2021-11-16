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
#occumb <- function(formula_phi = ~ 1,
#                   formula_theta = ~ 1,
#                   formula_psi = ~ 1,
#                   formula_phi_shared = NULL,
#                   formula_theta_shared = NULL,
#                   formula_psi_shared = NULL,
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
#    margs <- set_modargs(formula_phi,
#                         formula_theta,
#                         formula_psi,
#                         formula_phi_shared,
#                         formula_theta_shared,
#                         formula_psi_shared,
#                         data)

#    # Write model file
#    model <- tempfile()
#    writeLines(write_jags_model(margs$phi, margs$theta, margs$psi,
#                                margs$phi_shared,
#                                margs$phi_shared,
#                                margs$phi_shared), model)

    # Set data list
#    dat <- set_data(const, margs, prior_prec, prior_ulim)
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
#    inits <- set_inits(const, margs)
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
#    params <- set_params_monitored(margs$phi_shared,
#                                   margs$theta_shared,
#                                   margs$psi_shared)
    params <- c("Mu", "sigma", "rho", "alpha", "beta", "gamma",
                "phi", "theta", "psi", "z", "pi")

    # Run MCMC in JAGS
    fit <- jagsUI::jags(
#        dat, inits, params, model,
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
set_modargs <- function(formula_phi,
                        formula_theta,
                        formula_psi,
                        formula_phi_shared,
                        formula_theta_shared,
                        formula_psi_shared,
                        data) {

    phi_shared   <- !is.null(formula_phi_shared)
    theta_shared <- !is.null(formula_theta_shared)
    psi_shared   <- !is.null(formula_psi_shared)

    # theta_shared

    if (theta_shared) {
        # Stop when formula_theta_shared does not have an intercept
        check_intercept(formula_theta_shared, "theta_shared")
        # Stop when formula_theta_shared includes a term other than
        # site_cov, spec_cov, and their interaction
        check_wrong_terms(formula_theta_shared,
                          c(names(data@spec_cov),
                            names(data@site_cov),
                            names(data@repl_cov)),
                          "theta_shared")

        theta_shared_main_effects <- main_effects(terms(formula_theta_shared))

        # For theta = "ijk"
        if (any(theta_shared_main_effects %in% names(data@repl_cov))) {

        # For theta = "ij"
        } else if (any(theta_shared_main_effects %in% names(data@site_cov))) {

        # For theta = "i"
        } else {
            # Generate covariate objects
            for (n in seq_along(theta_shared_main_effects))
                eval(parse(text = sprintf("%s <- extract_covariate(theta_shared_main_effects[%s], data)", theta_shared_main_effects[n], n)))

            # Set design matrix
            dm <- set_design_matrix(formula_theta_shared, omit_intercept = TRUE)
            cov_theta_shared <- matrix(dm, nrow = dim(data@y)[1])
            colnames(cov_theta_shared) <- colnames(dm)
            M_theta_shared <- ncol(cov_theta_shared)
        }
    } else {
        cov_theta_shared <- NULL
        M_theta_shared   <- 0
    }



    # psi_shared

    if (psi_shared) {
        # Stop when formula_psi_shared does not have an intercept
        check_intercept(formula_psi_shared, "psi_shared")
        # Stop when formula_psi_shared includes a term other than
        # site_cov, spec_cov, and their interaction
        check_wrong_terms(formula_psi_shared,
                          c(names(data@spec_cov), names(data@site_cov)),
                          "psi_shared")

        psi_shared_main_effects <- main_effects(terms(formula_psi_shared))

        # For psi = "ij"
        if (any(psi_shared_main_effects %in% names(data@site_cov))) {

            # Generate covariate objects
            for (n in seq_along(psi_shared_main_effects)) {
                if (psi_shared_main_effects[n] %in% names(data@spec_cov))
                    eval(parse(text = sprintf("%s <- rep(extract_covariate(psi_shared_main_effects[%s], data), dim(data@y)[2])", psi_shared_main_effects[n], n)))
                if (psi_shared_main_effects[n] %in% names(data@site_cov))
                    eval(parse(text = sprintf("%s <- rep(extract_covariate(psi_shared_main_effects[%s], data), each = dim(data@y)[1])", psi_shared_main_effects[n], n)))
            }

            # Set design matrix
            dm <- set_design_matrix(formula_psi_shared, omit_intercept = TRUE)
            cov_psi_shared <- array(dm, c(dim(data@y)[1], dim(data@y)[2], ncol(dm)))
            dimnames(cov_psi_shared)[[3]] <- colnames(dm)
            M_psi_shared   <- dim(cov_psi_shared)[3]

        # For psi = "i"
        } else {
            # Generate covariate objects
            for (n in seq_along(psi_shared_main_effects))
                eval(parse(text = sprintf("%s <- extract_covariate(psi_shared_main_effects[%s], data)", psi_shared_main_effects[n], n)))

            # Set design matrix
            dm <- set_design_matrix(formula_psi_shared, omit_intercept = TRUE)
            cov_psi_shared <- matrix(dm, nrow = dim(data@y)[1])
            colnames(cov_psi_shared) <- colnames(dm)
            M_psi_shared   <- ncol(cov_psi_shared)
        }
    } else {
        cov_psi_shared <- NULL
        M_psi_shared   <- 0
    }

    M <- 0      # Order of species effects

    ## phi

    if (formula_phi == ~ 1) {
        cov_phi <- 1
        M       <- M + 1
        m_phi   <- 1
    } else {
        stop("Covariates in formula_phi are not yet supported, sorry.")
    }

    ## theta

    if (formula_theta == ~ 1) {
        cov_theta <- 1
        M         <- M + 1
        m_theta   <- seq(m_phi[length(m_phi)] + 1, m_phi[length(m_phi)] + 1)
    } else {
        # Stop when formula_theta does not have an intercept
        check_intercept(formula_theta, "theta")
        # Stop when formula_theta includes a term other than site/repl covs
        check_wrong_terms(formula_theta,
                          c(names(data@site_cov),
                            names(data@repl_cov)),
                          "theta")

        theta_main_effects <- main_effects(terms(formula_theta))

        # For theta = "ijk"
        if ((any(theta_main_effects %in% names(data@repl_cov)))) {
            # Generate covariate objects

            for (n in seq_along(theta_main_effects)) {
                if (theta_main_effects[n] %in% names(data@site_cov))
                    eval(parse(text = sprintf("%s <- rep(rep(extract_covariate(theta_main_effects[%s], data), each = dim(data@y)[1]), dim(data@y)[3])", theta_main_effects[n], n)))
                if (theta_main_effects[n] %in% names(data@repl_cov))
                    eval(parse(text = sprintf("%s <- rep(extract_covariate(theta_main_effects[%s], data), each = dim(data@y)[1])", theta_main_effects[n], n)))
            }

            # Set design matrix
            dm <- set_design_matrix(formula_theta)
            cov_theta <- array(dm, c(dim(data@y)[1], dim(data@y)[2], dim(data@y)[3], ncol(dm)))
            dimnames(cov_theta)[[4]] <- colnames(dm)
            M       <- M + dim(cov_theta)[4]
            m_theta <- seq(m_phi + 1, m_phi + dim(cov_theta)[4])

        # For theta = "ij"
        } else {
            # Generate covariate objects
            for (n in seq_along(theta_main_effects))
                eval(parse(text = sprintf("%s <- rep(extract_covariate(theta_main_effects[%s], data), each = dim(data@y)[1])", theta_main_effects[n], n)))

            # Set design matrix
            dm <- set_design_matrix(formula_theta)
            cov_theta <- array(dm, c(dim(data@y)[1], dim(data@y)[2], ncol(dm)))

            dimnames(cov_theta)[[3]] <- colnames(dm)
            M       <- M + dim(cov_theta)[3]
            m_theta <- seq(m_phi + 1, m_phi + dim(cov_theta)[3])
        }
    }

    ## psi

    if (formula_psi == ~ 1) {
        cov_psi <- 1
        M       <- M + 1
        m_psi   <- seq(m_theta[length(m_theta)] + 1, m_theta[length(m_theta)] + 1)

    # For psi = "ij"
    } else {
        # Stop when formula_psi does not have an intercept
        check_intercept(formula_psi, "psi")
        # Stop when formula_psi includes a term not in the site covariates
        check_wrong_terms(formula_psi, names(data@site_cov), "psi")

        # Generate covariate objects
        psi_main_effects <- main_effects(terms(formula_psi))
        for (n in seq_along(psi_main_effects))
            eval(parse(text = sprintf("%s <- rep(extract_covariate(psi_main_effects[%s], data), each = dim(data@y)[1])", psi_main_effects[n], n)))

        # Set design matrix
        dm <- set_design_matrix(formula_psi)
        cov_psi <- array(dm, c(dim(data@y)[1], dim(data@y)[2], ncol(dm)))
        dimnames(cov_psi)[[3]] <- colnames(dm)
        M       <- M + dim(cov_psi)[[3]]
        m_psi   <- seq(m_theta + 1, m_theta + dim(cov_psi)[[3]])
    }

    out <- list(phi              = "ijk",
                theta            = set_theta(formula_theta, formula_theta_shared, data),
                psi              = set_psi(formula_psi, formula_psi_shared, data),
                phi_shared       = phi_shared,
                theta_shared     = theta_shared,
                psi_shared       = psi_shared,
                M                = M,
                M_phi_shared     = 0,
                M_theta_shared   = M_theta_shared,
                M_psi_shared     = M_psi_shared,
                cov_phi          = cov_phi,
                cov_theta        = cov_theta,
                cov_psi          = cov_psi,
                cov_phi_shared   = 0,
                cov_theta_shared = cov_theta_shared,
                cov_psi_shared   = cov_psi_shared,
                m_phi            = m_phi,
                m_theta          = m_theta,
                m_psi            = m_psi)

    out
}
#set_modargs <- function(phi_formula, theta_formula, psi_formula, data) {
#    M    <- 0           # Order of species effects
#    code <- "000000"    # Model code
#
#    if (phi_formula == ~ 1) {
#        cov_phi <- 1
#        M       <- M + 1
#        m_phi   <- 1
#        substr(code, 1, 1) <- "1"
#    } else {
#        stop("Covariates in phi_formula are not yet supported, sorry.")
#    }
#
#    if (theta_formula == ~ 1) {
#        cov_theta <- 1
#        M         <- M + 1
#        m_theta   <- seq(m_phi + 1, m_phi + 1)
#        substr(code, 2, 2) <- "1"
#    } else {
#        stop("Covariates in theta_formula are not yet supported, sorry.")
#    }
#
#    if (psi_formula == ~ 1) {
#        cov_psi <- 1
#        M       <- M + 1
#        m_psi   <- seq(m_theta + 1, m_theta + 1)
#        substr(code, 3, 3) <- "1"
#    } else {
#        cov_psi <- stats::model.matrix(stats::as.formula(psi_formula),
#                                       data@site_cov)
#        M       <- M + ncol(cov_psi)
#        m_psi   <- seq(m_theta + 1, m_theta + ncol(cov_psi))
#        substr(code, 3, 3) <- "2"
#    }
#
#    out <- list(M         = M,
#                code      = code,
#                cov_phi   = cov_phi,
#                cov_theta = cov_theta,
#                cov_psi   = cov_psi,
#                m_phi     = m_phi,
#                m_theta   = m_theta,
#                m_psi     = m_psi)
#    out
#}

# Redefine the terms() function (!! DO NOT EXPORT !!)
terms <- function(formula) {
    if (is.null(formula))
        NULL
    else
        attr(stats::terms(formula), which = "term.labels")
}

main_effects <- function(terms) {
    if (is.null(terms))
        NULL
    else
        unique(unlist(strsplit(terms, split = ":")))
}

set_psi <- function(formula_psi, formula_psi_shared, data) {
    shared_main_effects <- main_effects(terms(formula_psi_shared))

    if (formula_psi == ~ 1) {
        if (is.null(shared_main_effects)) {
            psi <- "i"
        } else if (any(shared_main_effects %in% names(data@site_cov))) {
            psi <- "ij"
        } else {
            psi <- "i"
        }
    } else {
        psi <- "ij"
    }
}

set_theta <- function(formula_theta, formula_theta_shared, data) {
    m_eff  <- main_effects(terms(formula_theta))
    sm_eff <- main_effects(terms(formula_theta_shared))

    if (any(c(m_eff, sm_eff) %in% names(data@repl_cov))) {
        theta <- "ijk"
    } else {
        if (formula_theta == ~ 1) {
            if (is.null(sm_eff)) {
                theta <- "i"
            } else if (any(sm_eff %in% names(data@site_cov))) {
                theta <- "ij"
            } else {
                theta <- "i"
            }
        } else {
            theta <- "ij"
        }
    }
}

extract_covariate <- function(cov_name, data) {
    if (cov_name %in% names(data@spec_cov))
        cov_type <- "spec_cov"
    if (cov_name %in% names(data@site_cov))
        cov_type <- "site_cov"
    if (cov_name %in% names(data@repl_cov))
        cov_type <- "repl_cov"

    eval(parse(text = sprintf("data@%s$%s", cov_type, cov_name)))
}

`%!in%` <- Negate(`%in%`)

check_wrong_terms <- function(formula, correct_terms,
                              type = c("theta", "psi", "psi_shared", "theta_shared")) {
    test_terms <- terms(formula)
    wrong_terms <- main_effects(test_terms) %!in% correct_terms

    if (any(wrong_terms)) {
        if (type == "psi")
            stop(sprintf("Unexpected terms in formula_%s: %s
Note that only site covariates are allowed for formula_%s.",
                         type, test_terms[wrong_terms], type)) 
        if (type == "theta")
            stop(sprintf("Unexpected terms in formula_%s: %s
Note that species covariates are not allowed for formula_%s.",
                         type, test_terms[wrong_terms], type)) 
        if (type == "theta_shared")
            stop(sprintf("Unexpected terms in formula_%s: %s
Make sure they are found in either spec_cov, site_cov, or repl_cov.",
                         type, test_terms[wrong_terms])) 
        if (type == "psi_shared")
            stop(sprintf("Unexpected terms in formula_%s: %s
Note that only site covariates, species covariates, or their interactions are allowed for formula_%s.",
                         type, test_terms[wrong_terms], type)) 
    }
}

has_intercept <- function(formula) {
    as.logical(attributes(stats::terms(formula))$intercept)
}

check_intercept <- function(formula, type = c("psi", "psi_shared")) {
    if (!has_intercept(formula))
        stop(sprintf("No intercept in formula_%s: remove 0 or -1 from the formula",
                     type))
}

#xxxx <- function(formula, data, type = "psi_shared") {
#
#    # Generate covariate objects
#    m_eff <- main_effects(terms(formula))
#
#    if (type == "psi_shared") {
#        # For psi = "ij"
#        if (any(m_eff %in% names(data@site_cov))) {
#            for (n in seq_along(m_eff)) {
#                if (m_eff[n] %in% names(data@spec_cov))
#                    eval(parse(text = sprintf("%s <- rep(extract_covariate(m_eff[%s], data), dim(data@y)[2])", m_eff[n], n)))
#                if (m_eff[n] %in% names(data@site_cov))
#                    eval(parse(text = sprintf("%s <- rep(extract_covariate(m_eff[%s], data), each = dim(data@y)[1])", m_eff[n], n)))
#            }
#
#            tmp <- set_design_matrix(formula, type)
#            cov_psi_shared <- array(tmp, c(dim(data@y)[1], dim(data@y)[2], ncol(tmp)))
#            M_psi_shared   <- dim(cov_psi_shared)[3]
#
#        # For psi = "i"
#        } else {
#            for (n in seq_along(m_eff)) {
#                eval(parse(text = sprintf("%s <- extract_covariate(m_eff[%s], data)", m_eff[n], n)))
#            }
#
#            tmp <- set_design_matrix(formula, type)
#            cov_psi_shared <- matrix(tmp, nrow = dim(data@y)[1])
#            M_psi_shared   <- ncol(cov_psi_shared)
#        }
#    }
#
##    # Set design matrix
##    out <- stats::model.matrix(formula, parent.frame())
##    # When formula includes an intercept term, remove it and issue a warning
##    if (any(colnames(out) %in% "(Intercept)")) {
##        out <- out[, colnames(out) %!in% "(Intercept)"]
##        warning(sprintf("formula_%s should not include an intercept term: it will be removed.", type))
##    }
##
##    out
#}

set_design_matrix <- function(formula, omit_intercept = FALSE) {
    out <- stats::model.matrix(formula, parent.frame())

    if (omit_intercept) {
        if (sum(colnames(out) %!in% "(Intercept)") > 1) {
            out <- out[, colnames(out) %!in% "(Intercept)"]
        } else {
            cn <- colnames(out)[colnames(out) %!in% "(Intercept)"]
            out <- matrix(out[, colnames(out) %!in% "(Intercept)"], ncol = 1)
            colnames(out) <- cn
        }
    }

    out
}


# Auto-generate JAGS model code
write_jags_model <- function(phi, theta, psi, phi_shared, theta_shared, psi_shared) {

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

    if (phi_shared) {
        if (phi == "i")
            model <- c(model,
                       "        log(phi[i]) <- inprod(alpha[i, ], cov_phi[]) + inprod(alpha_shared[], cov_phi_shared[i, ])")
        if (phi == "ij")
            model <- c(model,
                       "        log(phi[i, j]) <- inprod(alpha[i, ], cov_phi[i, j, ]) + inprod(alpha_shared[], cov_phi_shared[i, j, ])")
        if (phi == "ijk")
            model <- c(model,
                       "        log(phi[i, j, k]) <- inprod(alpha[i, ], cov_phi[i, j, k, ]) + inprod(alpha_shared[], cov_phi_shared[i, j, k, ])")
    } else {
        if (phi == "i")
            model <- c(model,
                       "        log(phi[i]) <- inprod(alpha[i, ], cov_phi[])")
        if (phi == "ij")
            model <- c(model,
                       "        log(phi[i, j]) <- inprod(alpha[i, ], cov_phi[i, j, ])")
        if (phi == "ijk")
            model <- c(model,
                       "        log(phi[i, j, k]) <- inprod(alpha[i, ], cov_phi[i, j, k, ])")
    }

    if (theta_shared) {
        if (theta == "i")
            model <- c(model,
                       "        logit(theta[i]) <- inprod(beta[i, ], cov_theta[]) + inprod(beta_shared[], cov_theta_shared[i, ])")
        if (theta == "ij")
            model <- c(model,
                       "        logit(theta[i, j]) <- inprod(beta[i, ], cov_theta[i, j, ]) + inprod(beta_shared[], cov_theta_shared[i, j, ])")
        if (theta == "ijk")
            model <- c(model,
                       "        logit(theta[i, j, k]) <- inprod(beta[i, ], cov_theta[i, j, k, ]) + inprod(beta_shared[], cov_theta_shared[i, j, k, ])")
    } else {
        if (theta == "i")
            model <- c(model,
                       "        logit(theta[i]) <- inprod(beta[i, ], cov_theta[])")
        if (theta == "ij")
            model <- c(model,
                       "        logit(theta[i, j]) <- inprod(beta[i, ], cov_theta[i, j, ])")
        if (theta == "ijk")
            model <- c(model,
                       "        logit(theta[i, j, k]) <- inprod(beta[i, ], cov_theta[i, j, k, ])")
    }

    if (psi_shared) {
        if (psi == "i")
            model <- c(model,
                       "        logit(psi[i]) <- inprod(gamma[i, ], cov_psi[]) + inprod(gamma_shared[], cov_psi_shared[i, ])")
        if (psi == "ij")
            model <- c(model,
                       "        logit(psi[i, j]) <- inprod(gamma[i, ], cov_psi[i, j, ]) + inprod(gamma_shared[], cov_psi_shared[i, j, ])")
    } else {
        if (psi == "i")
            model <- c(model,
                       "        logit(psi[i]) <- inprod(gamma[i, ], cov_psi[])")
        if (psi == "ij")
            model <- c(model,
                       "        logit(psi[i, j]) <- inprod(gamma[i, ], cov_psi[i, j, ])")
    }

    model <- c(model,
               readLines(system.file("jags",
                                     "occumb_template5.jags",
                                     package = "occumb")))

    if (phi_shared)
        model <- c(model,
                   "    for (m in 1:M_phi_shared) {",
                   "        alpha_shared[m] ~ dnorm(0, prior_prec)",
                   "    }")
    if (theta_shared)
        model <- c(model,
                   "    for (m in 1:M_theta_shared) {",
                   "        beta_shared[m] ~ dnorm(0, prior_prec)",
                   "    }")
    if (psi_shared)
        model <- c(model,
                   "    for (m in 1:M_psi_shared) {",
                   "        gamma_shared[m] ~ dnorm(0, prior_prec)",
                   "    }")

    model <- c(model, "}", "")

    model
}

# Auto-generate the data list
set_data <- function(const, margs, prior_prec, prior_ulim) {
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

    if (margs$phi_shared)
        dat <- c(dat,
                 cov_phi_shared = margs$cov_phi_shared,
                 M_phi_shared   = margs$M_phi_shared)
    if (margs$theta_shared)
        dat <- c(dat,
                 cov_theta_shared = margs$cov_theta_shared,
                 M_theta_shared   = margs$M_theta_shared)
    if (margs$psi_shared)
        dat <- c(dat,
                 cov_psi_shared = margs$cov_psi_shared,
                 M_psi_shared   = margs$M_psi_shared)

    dat
}

# Auto-generate the initial value function
set_inits <- function(const, margs) {
    eval(parse(text = inits_code(const, margs)))
}

# Auto-generate codes for the initial value function
inits_code <- function(const, margs) {
    out <- c(
        "function() {",
        "    list(z = matrix(1, const$I, const$J),",
        "         u = array(1, dim = c(const$I, const$J, const$K)),",
        "         x = array(stats::rnorm(const$I * const$J * const$K,",
        "                                mean = 1, sd = 0.1),",
        "                   dim = c(const$I, const$J, const$K)),")

    if (margs$phi_shared)
        out <- c(out, "         alpha_shared = stats::rnorm(margs$M_phi_shared, sd = 0.1),")
    if (margs$theta_shared)
        out <- c(out, "         beta_shared  = stats::rnorm(margs$M_theta_shared, sd = 0.1),")
    if (margs$psi_shared)
        out <- c(out, "         gamma_shared = stats::rnorm(margs$M_psi_shared, sd = 0.1),")

    out <- c(out,
        "         spec_eff = matrix(stats::rnorm(const$I * margs$M, sd = 0.1),",
        "                           const$I, margs$M),",
        "         Mu       = stats::rnorm(margs$M, sd = 0.1),",
        "         sigma    = stats::rnorm(margs$M, mean = 1, sd = 0.1),",
        "         rho      = matrix(stats::rnorm(margs$M^2, sd = 0.1),",
        "                           margs$M, margs$M))",
        "}")

    out
}

# Auto-generate the list of parameters monitored
set_params_monitored <- function(phi_shared, theta_shared, psi_shared) {
    params <- c("Mu", "sigma", "rho", "alpha", "beta", "gamma")

    if (phi_shared)
        params <- c(params, "alpha_shared")
    if (theta_shared)
        params <- c(params, "beta_shared")
    if (psi_shared)
        params <- c(params, "gamma_shared")

    params <- c(params, c("phi", "theta", "psi", "z", "pi"))

    params
}

