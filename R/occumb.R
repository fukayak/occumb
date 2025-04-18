#' @include classes.R occumbData.R
NULL

# Class for model-fit results of occumb
setClass("occumbFit", slots = c(fit = "jagsUI",
                                data = "occumbData",
                                occumb_args = "list"))

#' @title Model-fitting function.
#' @description \code{occumb()} fits the multispecies site-occupancy model for eDNA
#'  metabarcoding (Fukaya et al. 2022) and returns a model-fit object containing
#'  posterior samples.
#' @details
#'  \code{occumb()} allows the fitting of a range of multispecies site
#'  occupancy models, including covariates at different levels of the data
#'  generation process.
#'  The most general form of the model can be written as follows (the notation
#'  follows that of the original article; see References).
#'
#'  Sequence read counts:
#'  \deqn{(y_{1jk}, ..., y_{Ijk}) \sim \textrm{Multinomial}((\pi_{1jk}, ...,  \pi_{Ijk}), N_{jk}),}{y[1:I, j, k] ~ Multinomial(pi[1:I, j, k], N[j, k]),}
#'  \deqn{\pi_{ijk} = \frac{u_{ijk}r_{ijk}}{\sum_m u_{mjk}r_{mjk}},}{pi[i, j, k] = (u[i, j, k] * r[i, j, k]) / sum(u[1:I, j, k] * r[1:I, j, k]),}
#'
#'  Relative frequency of species sequences:
#'  \deqn{r_{ijk} \sim \textrm{Gamma}(\phi_{ijk}, 1),}{r[i, j, k] ~ Gamma(phi[i, j, k], 1),}
#'
#'  Capture of species sequences:
#'  \deqn{u_{ijk} \sim \textrm{Bernoulli}(z_{ij}\theta_{ijk}),}{u[i, j, k] ~ Bernoulli(z[i, j] * theta[i, j, k]),}
#'
#'  Site occupancy of species:
#'  \deqn{z_{ij} \sim \textrm{Bernoulli}(\psi_{ij}),}{z[i, j] ~ Bernoulli(psi[i, j]),}
#'  where the variations of \eqn{\phi}{phi}, \eqn{\theta}{theta}, and
#'  \eqn{\psi}{psi} are modeled by specifying model formulas in
#'  \code{formula_phi}, \code{formula_theta}, \code{formula_psi},
#'  \code{formula_phi_shared}, \code{formula_theta_shared}, and
#'  \code{formula_psi_shared.}
#'  Each parameter may have species-specific effects and effects that are common
#'  across species, where the former is specified by \code{formula_phi},
#'  \code{formula_theta}, and \code{formula_psi}, whereas
#'  \code{formula_phi_shared}, \code{formula_theta_shared}, and
#'  \code{formula_psi_shared} specify the latter.
#'  As species-specific intercepts are specified by default, the intercept
#'  terms in \code{formula_phi_shared}, \code{formula_theta_shared}, and
#'  \code{formula_psi_shared} are always ignored.
#'  Covariate terms must be found in the names of the list elements stored
#'  in the \code{spec_cov}, \code{site_cov}, or \code{repl_cov} slots in the dataset
#'  object provided with the \code{data} argument.
#'  Covariates are modeled using the log link function for \eqn{\phi}{phi}
#'  and logit link function for \eqn{\theta}{theta} and \eqn{\psi.}{psi.}
#'
#'  The two arguments, \code{prior_prec} and \code{prior_ulim}, control the
#'  prior distribution of parameters. For the community-level average of
#'  species-specific effects and effects common across species, a normal
#'  prior distribution with a mean of 0 and precision (i.e., the inverse of the
#'  variance) \code{prior_prec} is specified. For the standard deviation of
#'  species-specific effects, a uniform prior distribution with a lower limit of
#'  zero and an upper limit of \code{prior_ulim} is specified. For the correlation
#'  coefficient of species-specific effects, a uniform prior distribution in the
#'  range of \eqn{-}1 to 1 is specified by default.
#'
#'  See \href{https://fukayak.github.io/occumb/articles/model_specification.html}{the package vignette}
#'  for details on the model specifications in \code{occumb()}.
#'
#'  The \code{data} argument requires a dataset object to be generated using
#'  \code{ocumbData()}; see the document of \code{\link{occumbData}()}.
#'
#'  The model is fit using the \code{\link[jagsUI]{jags}()} function of the
#'  \href{https://cran.r-project.org/package=jagsUI}{jagsUI} package, where
#'  Markov chain Monte Carlo methods are used to
#'  obtain posterior samples of the parameters and latent variables.
#'  Arguments \code{n.chains}, \code{n.adapt}, \code{n.burnin}, \code{n.thin},
#'  \code{n.iter}, and \code{parallel} are passed on to arguments of the
#'  same name in the \code{\link[jagsUI]{jags}()} function.
#'  See the document of \href{https://cran.r-project.org/package=jagsUI}{jagsUI}'s
#'  \code{\link[jagsUI]{jags}()} function for details.
#' @param formula_phi A right-hand side formula describing species-specific
#'        effects of sequence relative dominance (\eqn{\phi}).
#' @param formula_theta A right-hand side formula describing species-specific
#'        effects of sequence capture probability (\eqn{\theta}).
#' @param formula_psi A right-hand side formula describing species-specific
#'        effects of occupancy probability (\eqn{\psi}).
#' @param formula_phi_shared A right-hand side formula describing effects of
#'        sequence relative dominance (\eqn{\phi}) that are common across
#'        species. The intercept term is ignored (see Details).
#' @param formula_theta_shared A right-hand side formula describing effects of
#'        sequence capture probability (\eqn{\theta}) that are common across
#'        species.
#'        The intercept term is ignored (see Details).
#' @param formula_psi_shared A right-hand side formula describing effects of
#'        occupancy probability (\eqn{\psi}) that are common across species.
#'        The intercept term is ignored (see Details).
#' @param prior_prec Precision of normal prior distribution for the
#'        community-level average of species-specific parameters and effects
#'        common across species.
#' @param prior_ulim Upper limit of uniform prior distribution for the
#'        standard deviation of species-specific parameters.
#' @param data A dataset supplied as an \code{\link{occumbData}} class object.
#' @param n.chains Number of Markov chains to run.
#' @param n.adapt Number of iterations to run in the JAGS adaptive phase.
#' @param n.burnin Number of iterations at the beginning of the chain to discard.
#' @param n.thin Thinning rate. Must be a positive integer.
#' @param n.iter Total number of iterations per chain (including burn-in).
#' @param parallel If TRUE, run MCMC chains in parallel on multiple CPU cores.
#' @param ... Additional arguments passed to \code{\link[jagsUI]{jags}()} function.
#' @return  An S4 object of the \code{occumbFit} class containing the results of
#'          the model fitting and the supplied dataset.
#' @section References:
#'      K. Fukaya, N. I. Kondo, S. S. Matsuzaki and T. Kadoya (2022)
#'      Multispecies site occupancy modelling and study design for spatially
#'      replicated environmental DNA metabarcoding. \emph{Methods in Ecology
#'      and Evolution} \strong{13}:183--193.
#'      \doi{10.1111/2041-210X.13732}
#' @examples
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
#' \donttest{
#' # Fitting a null model (includes only species-specific intercepts)
#' res0 <- occumb(data = data)
#'
#' # Add species-specific effects of site covariates in occupancy probabilities
#' res1 <- occumb(formula_psi = ~ cov2, data = data)        # Continuous covariate
#' res2 <- occumb(formula_psi = ~ cov3, data = data)        # Categorical covariate
#' res3 <- occumb(formula_psi = ~ cov2 * cov3, data = data) # Interaction
#' # Add species covariate in the three parameters
#' # Note that species covariates are modeled as common effects
#' res4 <- occumb(formula_phi_shared = ~ cov1, data = data)   # phi
#' res5 <- occumb(formula_theta_shared = ~ cov1, data = data) # theta
#' res6 <- occumb(formula_psi_shared = ~ cov1, data = data)   # psi
#' # Add replicate covariates
#' # Note that replicate covariates can only be specified for theta and phi
#' res7 <- occumb(formula_phi = ~ cov4, data = data)   # phi
#' res8 <- occumb(formula_theta = ~ cov4, data = data) # theta
#' # Specify the prior distribution and MCMC settings explicitly
#' res9 <- occumb(data = data, prior_prec = 1E-2, prior_ulim = 1E2,
#'                n.chains = 1, n.burnin = 1000, n.thin = 1, n.iter = 2000)
#' res10 <- occumb(data = data, parallel = TRUE, n.cores = 2) # Run MCMC in parallel
#' }
#' @export
occumb <- function(formula_phi = ~ 1,
                   formula_theta = ~ 1,
                   formula_psi = ~ 1,
                   formula_phi_shared = ~ 1,
                   formula_theta_shared = ~ 1,
                   formula_psi_shared = ~ 1,
                   prior_prec = 1E-4,
                   prior_ulim = 1E4,
                   data,
                   n.chains = 4,
                   n.adapt = NULL,
                   n.burnin = 10000,
                   n.thin = 10,
                   n.iter = 20000,
                   parallel = FALSE,
                   ...) {

  # Validate arguments
  check_args_occumb(data, formula_phi, formula_theta, formula_psi,
                    formula_phi_shared, formula_theta_shared,
                    formula_psi_shared, prior_prec, prior_ulim)

  # Set constants
  const <- set_const(data)

  # Define arguments for the specified model
  margs <- set_modargs(formula_phi,
                       formula_theta,
                       formula_psi,
                       formula_phi_shared,
                       formula_theta_shared,
                       formula_psi_shared,
                       data)

  # Set data list
  dat <- set_data(const, margs, prior_prec, prior_ulim)

  # Set initial values
  inits <- set_inits(const, margs, dat)

  # Set parameters monitored
  params <- set_params_monitored(margs$phi_shared,
                                 margs$theta_shared,
                                 margs$psi_shared)

  # Write model file
  model <- tempfile()
  writeLines(write_jags_model(margs$phi, margs$theta, margs$psi,
                              margs$phi_shared,
                              margs$theta_shared,
                              margs$psi_shared), model)

  # Run MCMC in JAGS
  fit <- jagsUI::jags(dat, inits, params, model,
                      n.chains = n.chains,
                      n.adapt  = n.adapt,
                      n.iter   = n.iter,
                      n.burnin = n.burnin,
                      n.thin   = n.thin,
                      parallel = parallel, ...)

  # Output
  out <- methods::new(
    "occumbFit", fit = fit, data = data,
    occumb_args = list(
      formula_phi          = paste(as.character(formula_phi),
                                   collapse = " "),
      formula_theta        = paste(as.character(formula_theta),
                                   collapse = " "),
      formula_psi          = paste(as.character(formula_psi),
                                   collapse = " "),
      formula_phi_shared   = paste(as.character(formula_phi_shared),
                                   collapse = " "),
      formula_theta_shared = paste(as.character(formula_theta_shared),
                                   collapse = " "),
      formula_psi_shared   = paste(as.character(formula_psi_shared),
                                   collapse = " "),
      prior_prec           = prior_prec,
      prior_ulim           = prior_ulim
    )
  )
  out
}

# Validation for the inputs
check_args_occumb <- function(data,
                              formula_phi,
                              formula_theta,
                              formula_psi,
                              formula_phi_shared,
                              formula_theta_shared,
                              formula_psi_shared,
                              prior_prec,
                              prior_ulim) {
  # Check data
  if (!inherits(data, "occumbData"))
    stop("An occumbData class object is expected for data\n")

  # Check formulas
  formulas <- c("formula_phi",
                "formula_theta",
                "formula_psi",
                "formula_phi_shared",
                "formula_theta_shared",
                "formula_psi_shared")
  bad_formula <- rep(FALSE, 6)
  bad_formula[1] <- !inherits(formula_phi, "formula")
  bad_formula[2] <- !inherits(formula_theta, "formula")
  bad_formula[3] <- !inherits(formula_psi, "formula")
  bad_formula[4] <- !inherits(formula_phi_shared, "formula")
  bad_formula[5] <- !inherits(formula_theta_shared, "formula")
  bad_formula[6] <- !inherits(formula_psi_shared, "formula")
  if (any(bad_formula))
    stop(sprintf("Formula is expected for: %s\n",
                 paste(formulas[bad_formula], collapse = ", ")))

  check_intercept(formula_phi, "phi")
  check_intercept(formula_theta, "theta")
  check_intercept(formula_psi, "psi")
  check_intercept(formula_phi_shared, "phi_shared")
  check_intercept(formula_theta_shared, "theta_shared")
  check_intercept(formula_psi_shared, "psi_shared")

  # Check priors
  if (is.infinite(prior_prec))
    stop("'prior_prec' should have a finite value")
  if (!(prior_prec > 0))
    stop("'prior_prec' should have a positive value")
  if (!(prior_ulim > 0))
    stop("'prior_ulim' should have a positive value")
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

  check_wrong_terms(data, formula_phi, "phi")
  check_wrong_terms(data, formula_theta, "theta")
  check_wrong_terms(data, formula_psi, "psi")
  check_wrong_terms(data, formula_phi_shared, "phi_shared")
  check_wrong_terms(data, formula_theta_shared, "theta_shared")
  check_wrong_terms(data, formula_psi_shared, "psi_shared")

  phi_shared   <- formula_phi_shared   != ~ 1
  theta_shared <- formula_theta_shared != ~ 1
  psi_shared   <- formula_psi_shared   != ~ 1

  # phi
  list_covs_phi  <- set_covariates(data, formula_phi,
                                   formula_phi_shared, "phi")
  cov_phi        <- list_covs_phi$cov
  M_phi          <- list_covs_phi$M
  cov_phi_shared <- list_covs_phi$cov_shared
  M_phi_shared   <- list_covs_phi$M_shared
  type_phi       <- list_covs_phi$type

  # theta
  list_covs_theta  <- set_covariates(data, formula_theta,
                                     formula_theta_shared, "theta")
  cov_theta        <- list_covs_theta$cov
  M_theta          <- list_covs_theta$M
  cov_theta_shared <- list_covs_theta$cov_shared
  M_theta_shared   <- list_covs_theta$M_shared
  type_theta       <- list_covs_theta$type

  # psi
  list_covs_psi  <- set_covariates(data, formula_psi,
                                   formula_psi_shared, "psi")
  cov_psi        <- list_covs_psi$cov
  M_psi          <- list_covs_psi$M
  cov_psi_shared <- list_covs_psi$cov_shared
  M_psi_shared   <- list_covs_psi$M_shared
  type_psi       <- list_covs_psi$type

  M <- M_phi + M_theta + M_psi    # Order of species effects
  m_phi   <- seq(1, M_phi)
  m_theta <- seq(m_phi[length(m_phi)] + 1, m_phi[length(m_phi)] + M_theta)
  m_psi   <- seq(m_theta[length(m_theta)] + 1, m_theta[length(m_theta)] + M_psi)

  out <- list(phi          = type_phi,
              theta        = type_theta,
              psi          = type_psi,
              phi_shared   = phi_shared,
              theta_shared = theta_shared,
              psi_shared   = psi_shared,
              M            = M,
              cov_phi      = cov_phi,
              cov_theta    = cov_theta,
              cov_psi      = cov_psi,
              m_phi        = m_phi,
              m_theta      = m_theta,
              m_psi        = m_psi)

  if (phi_shared) {
    out$M_phi_shared   <- M_phi_shared
    out$cov_phi_shared <- cov_phi_shared
  }
  if (theta_shared) {
    out$M_theta_shared   <- M_theta_shared
    out$cov_theta_shared <- cov_theta_shared
  }
  if (psi_shared) {
    out$M_psi_shared   <- M_psi_shared
    out$cov_psi_shared <- cov_psi_shared
  }

  out
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

  if (margs$phi_shared) {
    dat$cov_phi_shared <- margs$cov_phi_shared
    dat$M_phi_shared   <- margs$M_phi_shared
  }
  if (margs$theta_shared) {
    dat$cov_theta_shared <- margs$cov_theta_shared
    dat$M_theta_shared   <- margs$M_theta_shared
  }
  if (margs$psi_shared) {
    dat$cov_psi_shared <- margs$cov_psi_shared
    dat$M_psi_shared   <- margs$M_psi_shared
  }

  dat
}

# Auto-generate the initial value function
set_inits <- function(const, margs, dat) {
  eval(parse(text = inits_code(const, margs, dat)))
}

# Auto-generate codes for the initial value function
inits_code <- function(const, margs, dat) {
  out <- c(
           "function() {",
           "    list(z = matrix(1, const$I, const$J),",
           "         u = array(1, dim = c(const$I, const$J, const$K)),",
           "         r = array(stats::rnorm(const$I * const$J * const$K,",
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
           "         Mu       = stats::rnorm(margs$M, sd = 0.1),")

  if (dat$prior_ulim > 10) {
    out <- c(out,
             "         sigma    = stats::rnorm(margs$M, mean = 1, sd = 0.1))")
  } else {
    out <- c(out,
             sprintf("         sigma    = stats::runif(margs$M, min = 0, max = %s))",
                     dat$prior_ulim))
  }

  out <- c(out, "}")

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

# Auto-generate JAGS model code
write_jags_model <- function(phi, theta, psi,
                             phi_shared, theta_shared, psi_shared) {

  model <- readLines(system.file("jags",
                                 "occumb_template1.jags",
                                 package = "occumb"))

  if (phi == "i")
    model <- c(model,
               "                r[i, j, k] ~ dgamma(phi[i], 1)")
  if (phi == "ij")
    model <- c(model,
               "                r[i, j, k] ~ dgamma(phi[i, j], 1)")
  if (phi == "ijk")
    model <- c(model,
               "                r[i, j, k] ~ dgamma(phi[i, j, k], 1)")

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
                 "        for (j in 1:J) {",
                 "            log(phi[i, j]) <- inprod(alpha[i, ], cov_phi[j, ]) + inprod(alpha_shared[], cov_phi_shared[i, j, ])",
                 "        }")
    if (phi == "ijk")
      model <- c(model,
                 "        for (j in 1:J) {",
                 "            for (k in 1:K) {",
                 "                log(phi[i, j, k]) <- inprod(alpha[i, ], cov_phi[j, k, ]) + inprod(alpha_shared[], cov_phi_shared[i, j, k, ])",
                 "            }",
                 "        }")
  } else {
    if (phi == "i")
      model <- c(model,
                 "        log(phi[i]) <- inprod(alpha[i, ], cov_phi[])")
    if (phi == "ij")
      model <- c(model,
                 "        for (j in 1:J) {",
                 "            log(phi[i, j]) <- inprod(alpha[i, ], cov_phi[j, ])",
                 "        }")
    if (phi == "ijk")
      model <- c(model,
                 "        for (j in 1:J) {",
                 "            for (k in 1:K) {",
                 "                log(phi[i, j, k]) <- inprod(alpha[i, ], cov_phi[j, k, ])",
                 "            }",
                 "        }")
  }

  if (theta_shared) {
    if (theta == "i")
      model <- c(model,
                 "        logit(theta[i]) <- inprod(beta[i, ], cov_theta[]) + inprod(beta_shared[], cov_theta_shared[i, ])")
    if (theta == "ij")
      model <- c(model,
                 "        for (j in 1:J) {",
                 "            logit(theta[i, j]) <- inprod(beta[i, ], cov_theta[j, ]) + inprod(beta_shared[], cov_theta_shared[i, j, ])",
                 "        }")
    if (theta == "ijk")
      model <- c(model,
                 "        for (j in 1:J) {",
                 "            for (k in 1:K) {",
                 "                logit(theta[i, j, k]) <- inprod(beta[i, ], cov_theta[j, k, ]) + inprod(beta_shared[], cov_theta_shared[i, j, k, ])",
                 "            }",
                 "        }")
  } else {
    if (theta == "i")
      model <- c(model,
                 "        logit(theta[i]) <- inprod(beta[i, ], cov_theta[])")
    if (theta == "ij")
      model <- c(model,
                 "        for (j in 1:J) {",
                 "            logit(theta[i, j]) <- inprod(beta[i, ], cov_theta[j, ])",
                 "        }")
    if (theta == "ijk")
      model <- c(model,
                 "        for (j in 1:J) {",
                 "            for (k in 1:K) {",
                 "                logit(theta[i, j, k]) <- inprod(beta[i, ], cov_theta[j, k, ])",
                 "            }",
                 "        }")
  }

  if (psi_shared) {
    if (psi == "i")
      model <- c(model,
                 "        logit(psi[i]) <- inprod(gamma[i, ], cov_psi[]) + inprod(gamma_shared[], cov_psi_shared[i, ])")
    if (psi == "ij")
      model <- c(model,
                 "        for (j in 1:J) {",
                 "            logit(psi[i, j]) <- inprod(gamma[i, ], cov_psi[j, ]) + inprod(gamma_shared[], cov_psi_shared[i, j, ])",
                 "        }")
  } else {
    if (psi == "i")
      model <- c(model,
                 "        logit(psi[i]) <- inprod(gamma[i, ], cov_psi[])")
    if (psi == "ij")
      model <- c(model,
                 "        for (j in 1:J) {",
                 "            logit(psi[i, j]) <- inprod(gamma[i, ], cov_psi[j, ])",
                 "        }")
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


# Auxiliary functions for set_modargs() ----------------------------------------
set_psi <- function(formula_psi, formula_psi_shared, data) {
  shared_main_effects <- all.vars(formula_psi_shared)

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

  psi
}

set_phi_theta <- function(formula, formula_shared, data) {
  m_eff  <- all.vars(formula)
  sm_eff <- all.vars(formula_shared)

  if (any(c(m_eff, sm_eff) %in% names(data@repl_cov))) {
    out <- "ijk"
  } else {
    if (formula == ~ 1) {
      if (is.null(sm_eff)) {
        out <- "i"
      } else if (any(sm_eff %in% names(data@site_cov))) {
        out <- "ij"
      } else {
        out <- "i"
      }
    } else {
      out <- "ij"
    }
  }

  out
}

check_intercept <- function(formula, type) {
  if ("." %in% all.vars(formula)) return()
  if (!has_intercept(formula))
    stop(sprintf("No intercept in formula_%s: remove 0 or -1 from the formula\n",
                 type))
}

has_intercept <- function(formula) {
  as.logical(attributes(stats::terms(formula))$intercept)
}

check_wrong_terms <- function(data, formula, type) {

  type <- match.arg(type, c("phi", "phi_shared",
                            "theta", "theta_shared",
                            "psi", "psi_shared"))

  vars <- all.vars(formula)

  if ("." %in% vars)
    stop(sprintf("'.' cannot be used to specify covariates (formula_%s)", type))

  correct_terms <- if (type == "phi" || type == "theta") {
    c(names(data@site_cov), names(data@repl_cov))
  } else if (type == "psi") {
    names(data@site_cov)
  } else if (type == "phi_shared" || type == "theta_shared") {
    c(names(data@spec_cov), names(data@site_cov), names(data@repl_cov))
  } else if (type == "psi_shared") {
    c(names(data@spec_cov), names(data@site_cov))
  }

  wrong_terms <- vars %!in% correct_terms

  if (any(wrong_terms)) {
    if (type == "phi")
      stop(sprintf("Unexpected terms in formula_%s: %s
Make sure they are included in either 'site_cov' or 'repl_cov'.

For details on covariate specification in occumb(), see: https://fukayak.github.io/occumb/articles/model_specification.html#covariate-modeling-of-psi-theta-and-phi\n",
                   type,
                   knitr::combine_words(vars[wrong_terms], and = "")))
    if (type == "theta")
      stop(sprintf("Unexpected terms in formula_%s: %s
Make sure they are included in either 'site_cov' or 'repl_cov'.

For details on covariate specification in occumb(), see: https://fukayak.github.io/occumb/articles/model_specification.html#covariate-modeling-of-psi-theta-and-phi\n",
                   type,
                   knitr::combine_words(vars[wrong_terms], and = "")))
    if (type == "psi")
      stop(sprintf("Unexpected terms in formula_%s: %s
Make sure they are included in 'site_cov'.

For details on covariate specification in occumb(), see: https://fukayak.github.io/occumb/articles/model_specification.html#covariate-modeling-of-psi-theta-and-phi\n",
                   type,
                   knitr::combine_words(vars[wrong_terms], and = "")))
    if (type == "phi_shared")
      stop(sprintf("Unexpected terms in formula_%s: %s
Make sure they are included in either 'spec_cov', 'site_cov', or 'repl_cov'.

For details on covariate specification in occumb(), see: https://fukayak.github.io/occumb/articles/model_specification.html#covariate-modeling-of-psi-theta-and-phi\n",
                   type,
                   knitr::combine_words(vars[wrong_terms], and = "")))
    if (type == "theta_shared")
      stop(sprintf("Unexpected terms in formula_%s: %s
Make sure they are included in either 'spec_cov', 'site_cov', or 'repl_cov'.

For details on covariate specification in occumb(), see: https://fukayak.github.io/occumb/articles/model_specification.html#covariate-modeling-of-psi-theta-and-phi\n",
                   type,
                   knitr::combine_words(vars[wrong_terms], and = "")))
    if (type == "psi_shared")
      stop(sprintf("Unexpected terms in formula_%s: %s
Make sure they are included in either 'spec_cov' or 'site_cov'.

For details on covariate specification in occumb(), see: https://fukayak.github.io/occumb/articles/model_specification.html#covariate-modeling-of-psi-theta-and-phi\n",
                   type,
                   knitr::combine_words(vars[wrong_terms], and = "")))
  }
}

set_design_matrix <- function(formula, list_cov, omit_intercept = FALSE) {

  if (formula == ~ 1) {
    return(invisible())
  } else {
    out <- tryCatch({
      stats::model.matrix(formula, data = list_cov)
    }, warning = function(w) {
      formula_terms <- formula |>
        stats::terms() |>
        attr("variable")

      invalid <- formula_terms |>
        sapply(function(x) {
          inherits(x, "call") &&
            x[[1]] == "poly" &&
            eval(parse(text = sprintf("is.character(list_cov$%s)", all.vars(x))))
        }
        )

      if (any(invalid)) {
        stop(
          sprintf(
            "Numerical operations on character covariates '%s': make sure to avoid using operations like `poly` on characters\n",
            sapply(formula_terms, all.vars)[invalid]
          ),
          call. = FALSE
        )
      } else {
        stop(w, "\n", call. = FALSE)
      }
    })
  }

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

set_covariates <- function(data, formula, formula_shared, param) {

  param <- match.arg(param, c("phi", "theta", "psi"))

  ## phi and theta
  if (param != "psi" && formula_shared == ~ 1) {
    # For cases when formula_shared is not specified
    cov_shared <- M_shared <- NULL

    if (formula == ~ 1) {
      cov <- 1
      M <- 1
      type <- "i"
    } else {

      type <- set_phi_theta(formula, formula_shared, data)

      if (type == "ij") {
        dm <- set_design_matrix(formula, get_list_cov(data, "ij"))
        cov <- matrix(dm, nrow = dim(data@y)[2], ncol = ncol(dm))
        colnames(cov) <- colnames(dm)
        M <- ncol(cov)
      }

      if (type == "ijk") {
        dm <- set_design_matrix(formula, get_list_cov(data, "ijk"))
        cov <- array(dm, c(dim(data@y)[2], dim(data@y)[3], ncol(dm)))
        dimnames(cov)[[3]] <- colnames(dm)
        M <- dim(cov)[3]
      }
    }
  }
  if (param != "psi" && formula_shared != ~ 1) {
    # For cases when formula_shared is specified

    type <- set_phi_theta(formula, formula_shared, data)

    if (type == "i") {
      dm_shared <-
        set_design_matrix(formula_shared,
                          get_list_cov(data, "i", shared = TRUE),
                          omit_intercept = TRUE)
      cov <- 1
      M   <- 1
      cov_shared <- matrix(dm_shared, nrow = dim(data@y)[1])
      colnames(cov_shared) <- colnames(dm_shared)
      M_shared <- ncol(cov_shared)
    }

    if (type == "ij") {
      dm <- set_design_matrix(formula, get_list_cov(data, "ij"))
      dm_shared <-
        set_design_matrix(formula_shared,
                          get_list_cov(data, "ij", shared = TRUE),
                          omit_intercept = TRUE)
      if (formula == ~ 1) {
        cov <- matrix(1, nrow = dim(data@y)[2], ncol = 1)
        M   <- 1
      } else {
        cov <- matrix(dm, nrow = dim(data@y)[2], ncol = ncol(dm))
        colnames(cov) <- colnames(dm)
        M <- ncol(cov)
      }
      cov_shared <- array(dm_shared, c(dim(data@y)[1],
                                       dim(data@y)[2],
                                       ncol(dm_shared)))
      dimnames(cov_shared)[[3]] <- colnames(dm_shared)
      M_shared <- dim(cov_shared)[3]
    }

    if (type == "ijk") {
      dm <- set_design_matrix(formula, get_list_cov(data, "ijk"))
      dm_shared <-
        set_design_matrix(formula_shared,
                          get_list_cov(data, "ijk", shared = TRUE),
                          omit_intercept = TRUE)
      if (formula == ~ 1) {
        cov <- array(1, dim = c(dim(data@y)[2], dim(data@y)[3], 1))
        M   <- 1
      } else {
        cov <- array(dm,
                     c(dim(data@y)[2], dim(data@y)[3], ncol(dm)))
        dimnames(cov)[[3]] <- colnames(dm)
        M <- dim(cov)[3]
      }
      cov_shared <- array(dm_shared, c(dim(data@y)[1],
                                       dim(data@y)[2],
                                       dim(data@y)[3],
                                       ncol(dm_shared)))
      dimnames(cov_shared)[[4]] <- colnames(dm_shared)
      M_shared <- dim(cov_shared)[4]
    }
  }

  if (param == "psi" && formula_shared == ~ 1) {
    # For cases when formula_shared is not specified
    cov_shared <- M_shared <- NULL

    if (formula == ~ 1) {
      cov <- 1
      M <- 1
      type <- "i"
    } else {
      dm <- set_design_matrix(formula, get_list_cov(data, "ij"))
      cov <- matrix(dm, nrow =  dim(data@y)[2], ncol = ncol(dm))
      colnames(cov) <- colnames(dm)
      M <- ncol(cov)
      type <- "ij"
    }
  }
  if (param == "psi" && formula_shared != ~ 1) {
    # For cases when formula_shared is specified

    type <- set_psi(formula, formula_shared, data)

    if (type == "i") {
      dm_shared <-
        set_design_matrix(formula_shared,
                          get_list_cov(data, "i", shared = TRUE),
                          omit_intercept = TRUE)
      cov <- 1
      M   <- 1
      cov_shared <- matrix(dm_shared, nrow = dim(data@y)[1])
      colnames(cov_shared) <- colnames(dm_shared)
      M_shared <- ncol(cov_shared)
    }

    if (type == "ij") {
      dm <- set_design_matrix(formula, get_list_cov(data, "ij"))
      dm_shared <-
        set_design_matrix(formula_shared,
                          get_list_cov(data, "ij", shared = TRUE),
                          omit_intercept = TRUE)
      if (formula == ~ 1) {
        cov <- matrix(1, nrow = dim(data@y)[2], ncol = 1)
        M   <- 1
      } else {
        cov <- matrix(dm, nrow = dim(data@y)[2], ncol = ncol(dm))
        colnames(cov) <- colnames(dm)
        M <- ncol(cov)
      }
      cov_shared <- array(dm_shared, c(dim(data@y)[1],
                                       dim(data@y)[2],
                                       ncol(dm_shared)))
      dimnames(cov_shared)[[3]] <- colnames(dm_shared)
      M_shared <- dim(cov_shared)[3]
    }
  }

  out <- list(cov        = cov,
              M          = M,
              cov_shared = cov_shared,
              M_shared   = M_shared,
              type       = type)
  return(out)
}

get_list_cov <- function(data, type, shared = FALSE) {
  type <- match.arg(type, c("i", "ij", "ijk"))
  out <- list()

  if (!shared && type == "ij") {
    for (n in seq_along(names(data@site_cov))) {
      eval(parse(text = sprintf("out$%s <- data@site_cov[[%s]]",
                                names(data@site_cov)[n], n)))
    }
  }

  if (!shared && type == "ijk") {
    for (n in seq_along(names(data@site_cov))) {
      eval(parse(text = sprintf("out$%s <- rep(data@site_cov[[%s]], dim(data@y)[3])",
                                names(data@site_cov)[n], n)))
    }
    for (n in seq_along(names(data@repl_cov))) {
      eval(parse(text = sprintf("out$%s <- c(data@repl_cov[[%s]])",
                                names(data@repl_cov)[n], n)))
    }
  }

  if (shared && type == "i") {
    for (n in seq_along(names(data@spec_cov))) {
      eval(parse(text = sprintf("out$%s <- data@spec_cov[[%s]]",
                                names(data@spec_cov)[n], n)))
    }
  }

  if (shared && type == "ij") {
    for (n in seq_along(names(data@spec_cov))) {
      eval(parse(text = sprintf("out$%s <- rep(data@spec_cov[[%s]], dim(data@y)[2])",
                                names(data@spec_cov)[n], n)))
    }
    for (n in seq_along(names(data@site_cov))) {
      eval(parse(text = sprintf("out$%s <- rep(data@site_cov[[%s]], each = dim(data@y)[1])",
                                names(data@site_cov)[n], n)))
    }
  }

  if (shared && type == "ijk") {
    for (n in seq_along(names(data@spec_cov))) {
      eval(parse(text = sprintf("out$%s <- rep(rep(data@spec_cov[[%s]], dim(data@y)[2]), dim(data@y)[3])",
                                names(data@spec_cov)[n], n)))
    }
    for (n in seq_along(names(data@site_cov))) {
      eval(parse(text = sprintf("out$%s <- rep(rep(data@site_cov[[%s]], each = dim(data@y)[1]), dim(data@y)[3])",
                                names(data@site_cov)[n], n)))
    }
    for (n in seq_along(names(data@repl_cov))) {
      eval(parse(text = sprintf("out$%s <- rep(data@repl_cov[[%s]], each = dim(data@y)[1])",
                                names(data@repl_cov)[n], n)))
    }
  }

  return(out)
}


# -----------------------------------------------------------------------------


# Getter ----------------------------------------------------------------------
get_data <- function(occumbFit, variable) {
  eval(parse(text = paste0("occumbFit@data@", variable)))
}
get_modargs <- function(occumbFit) {
  set_modargs(stats::as.formula(occumbFit@occumb_args$formula_phi),
              stats::as.formula(occumbFit@occumb_args$formula_theta),
              stats::as.formula(occumbFit@occumb_args$formula_psi),
              stats::as.formula(occumbFit@occumb_args$formula_phi_shared),
              stats::as.formula(occumbFit@occumb_args$formula_theta_shared),
              stats::as.formula(occumbFit@occumb_args$formula_psi_shared),
              occumbFit@data)
}
# -----------------------------------------------------------------------------
