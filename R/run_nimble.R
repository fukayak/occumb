#' Run MCMC with the NIMBLE backend (internal)
#'
#' Internal worker used to fit the model with the NIMBLE engine. This function
#' builds a NIMBLE model, runs MCMC (optionally in parallel), summarizes the
#' posterior samples, and returns an object compatible with the
#' \pkg{jagsUI}-style fit object.
#'
#' @param data A named list containing the observed data and covariates for NIMBLE.
#' @param const A named list of model constants for NIMBLE.
#' @param inits A function that returns a named list of initial values for the
#'   model parameters (one set of initial values will be generated per chain).
#' @param params Character vector of node names to monitor.
#' @param model_code_strings Character vector of NIMBLE model code lines to be
#'   converted to \code{nimbleCode}.
#' @param model_file A string identifying the model file (kept for compatibility
#'   with other engines); stored in the returned object.
#' @param n.chains Number of MCMC chains.
#' @param n.iter Total number of MCMC iterations per chain.
#' @param n.burnin Number of burn-in iterations.
#' @param n.thin Thinning interval.
#' @param parallel Logical; if \code{TRUE}, run chains in parallel using the
#'   \pkg{parallel} package.
#' @param ... Additional control arguments:
#'   \describe{
#'     \item{\code{seed}}{See \code{nimble::runMCMC(setSeed = ...)}.
#'       \code{FALSE}: no seeding; \code{TRUE}: seed chain i with i;
#'       numeric vector (\code{length = n.chains}): per-chain seeds.}
#'     \item{\code{n.cores}}{Number of worker processes when \code{parallel=TRUE}.
#'       Defaults to \code{parallel::detectCores() - 1} (minimum 1), capped at
#'       \code{n.chains}.}
#'     \item{\code{store.data}}{Logical; if \code{TRUE}, store the input
#'       \code{data} and generated initial values in the returned object.}
#'     \item{\code{verbose}}{Logical; if not \code{NULL}, temporarily set
#'       \code{nimble} options \code{verbose} and \code{MCMCprogressBar}
#'       accordingly during execution.}
#'   }
#'
#' @return \pkg{jagsUI}-style fit object.
#' @keywords internal
run_nimble <- function(data, const, inits, params, model_code_strings, model_file,
                       n.chains, n.iter, n.burnin, n.thin, parallel, ...) {
  attach_nimble_package()

  # Check arguments
  model_code <- to_model_code(model_code_strings)
  dots_arguments <- list(...)
  n.cores <- dots_arguments$n.cores
  stopifnot(is.null(n.cores) || (is.numeric(n.cores) && length(n.cores) == 1L && n.cores > 0))
  seed <- dots_arguments$seed
  stopifnot(is.null(seed) || is.logical(seed) || is.numeric(seed))
  store.data <- dots_arguments$store.data
  stopifnot(is.null(store.data) || is.logical(store.data))
  if (is.null(store.data)) store.data <- FALSE
  verbose <- dots_arguments$verbose
  stopifnot(is.null(verbose) || is.logical(verbose))
  if (!is.null(verbose)) {
    verbose_old <- nimble::getNimbleOption("verbose")
    progress_bar_old <- nimble::getNimbleOption("MCMCprogressBar")
    nimble::nimbleOptions(verbose = verbose, MCMCprogressBar = verbose)
    on.exit({
      nimble::nimbleOptions(verbose = verbose_old)
      nimble::nimbleOptions(MCMCprogressBar = progress_bar_old)
    }, add = TRUE)
  }

  # Set constants
  const_nimble <- set_const_nimble(const, data)

  # Set data list
  data_nimble <- set_data_nimble(data)

  # Set initial values
  n_rho <- max(const_nimble$rho_index)
  inits_nimble <- set_inits_nimble(inits, seed, n.chains, n_rho)

  # Run MCMC in NIMBLE
  start_time <- Sys.time()
  if (parallel) {
    fit <- run_nimble_parallel(inits = inits_nimble, code = model_code,
                               const = const_nimble, data = data_nimble,
                               monitors = params,
                               n.iter, n.burnin, n.thin, n.chains,
                               n.cores = n.cores)
  } else {
    fit <- run_nimble_model(inits = inits_nimble, code = model_code,
                            const = const_nimble, data = data_nimble,
                            monitors = params,
                            n.iter, n.burnin, n.thin, n.chains)
  }
  elapsed_mins <- round(as.numeric(Sys.time() - start_time, units = "mins"),
                        digits = 3)
  nimble::messageIfVerbose("Summarizing MCMC samples...")
  fit <- nimbleSummary(fit)
  fit <- make_jagsui_compatible(fit)
  nimble::messageIfVerbose("Finished")

  fit
}

run_nimble_parallel <- function(inits, code, const, data, monitors,
                                n.iter, n.burnin, n.thin, n.chains, n.cores) {
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required. Please install it.", call. = FALSE)
  }
  if (is.null(n.cores)) {
    n.cores <- detect_cores_omit_one()
    n.cores <- min(n.cores, n.chains)
    nimble::messageIfVerbose(sprintf("[Note] Automatically setting n.cores = %d.", n.cores))
  }
  if (n.chains < n.cores) {
    n.cores <- n.chains
    nimble::messageIfVerbose("[Note] 'n.cores' exceeds 'n.chains'; reducing 'n.cores' to 'n.chains'.")
  }

  cluster <- parallel::makeCluster(n.cores)
  on.exit({
    tryCatch(parallel::stopCluster(cl = cluster), error = function(e) NULL)
    for (nd in cluster)
      tryCatch(close(nd$con), error = function(e) NULL)
  }, add = TRUE)

  parallel::clusterEvalQ(cluster, library(nimble))
  results <- parallel::parLapply(cl = cluster, X = inits,
                                 fun = run_nimble_model, code = code,
                                 const = const, data = data,
                                 monitors = monitors,
                                 n.iter = n.iter, n.burnin = n.burnin,
                                 n.thin = n.thin, n.chains = 1L,
                                 parallel = TRUE)
  parallel::stopCluster(cl = cluster)

  results
}

run_nimble_model <- function(inits, code, const, data, monitors,
                             n.iter, n.burnin, n.thin, n.chains,
                             parallel = FALSE) {
  if (parallel) {
    seed  <- inits$.RNG.seed
    inits <- inits[ls(inits)]
    inits_model <- inits
  } else {
    seed  <- vapply(inits, function(x) x$.RNG.seed, FUN.VALUE = numeric(1L))
    inits <- lapply(inits, function(x) x[ls(x)])
    inits_model <- inits[[1]]
  }

  model  <- nimble::nimbleModel(code = code, constants = const, data = data,
                                inits = inits_model, buildDerivs = TRUE)
  Cmodel <- nimble::compileNimble(model)
  conf   <- nimble::configureMCMC(Cmodel, monitors = monitors, print = FALSE)
  conf$replaceSamplers(target = "Mu", type = "barker", silent = TRUE)
  conf$replaceSamplers(target = "sigma", type = "barker", silent = TRUE)
  inds_r <- conf$findSamplersOnNodes("r")
  conf$removeSamplers(ind = inds_r)
  for (j in seq_len(const$J)) {
    for (k in seq_len(const$K)) {
      target <- sprintf("r[, %d, %d]", j, k)
      conf$addSampler(target = target, type = "barker", silent = TRUE)
    }
  }
  MCMC   <- nimble::buildMCMC(conf)
  CMCMC  <- nimble::compileNimble(MCMC, project = Cmodel)
  result <- nimble::runMCMC(CMCMC, niter = n.iter, nburnin = n.burnin,
                            thin = n.thin, nchains = n.chains, inits = inits,
                            setSeed = seed)
  result
}

# Auto-generate JAGS model code
write_nimble_model <- function(phi, theta, psi,
                               phi_shared, theta_shared, psi_shared,
                               M_cov_phi, M_cov_phi_shared,
                               M_cov_theta, M_cov_theta_shared,
                               M_cov_psi, M_cov_psi_shared) {

  model <- readLines(system.file("nimble",
                                 "occumb_template1.nimble",
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
             readLines(system.file("nimble",
                                   "occumb_template2.nimble",
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
    if (phi == "i") {
      if (M_cov_phi == 1) {
        term1 <- "alpha[i, 1] * cov_phi[1]"
      } else {
        term1 <- "inprod(alpha[i, 1:M_phi], cov_phi[1:M_phi])"
      }
      if (M_cov_phi_shared == 1) {
        term2 <- "alpha_shared[1] * cov_phi_shared[i, 1]"
      } else {
        term2 <- "inprod(alpha_shared[1:M_phi_shared], cov_phi_shared[i, 1:M_phi_shared])"
      }
      model <- c(model, paste0(
        "        log(phi[i]) <- ", term1, " + ", term2))
    } else if (phi == "ij") {
      if (M_cov_phi == 1) {
        term1 <- "alpha[i, 1] * cov_phi[j, 1]"
      } else {
        term1 <- "inprod(alpha[i, 1:M_phi], cov_phi[j, 1:M_phi])"
      }
      if (M_cov_phi_shared == 1) {
        term2 <- "alpha_shared[1] * cov_phi_shared[i, j, 1]"
      } else {
        term2 <- "inprod(alpha_shared[1:M_phi_shared], cov_phi_shared[i, j, 1:M_phi_shared])"
      }
      model <- c(model,
                 "        for (j in 1:J) {", paste0(
                   "            log(phi[i, j]) <- ", term1, " + ", term2),
                 "        }")
    } else if (phi == "ijk") {
      if (M_cov_phi == 1) {
        term1 <- "alpha[i, 1] * cov_phi[j, k, 1]"
      } else {
        term1 <- "inprod(alpha[i, 1:M_phi], cov_phi[j, k, 1:M_phi])"
      }
      if (M_cov_phi_shared == 1) {
        term2 <- "alpha_shared[1] * cov_phi_shared[i, j, k, 1]"
      } else {
        term2 <- "inprod(alpha_shared[1:M_phi_shared], cov_phi_shared[i, j, k, 1:M_phi_shared])"
      }
      model <- c(model,
                 "        for (j in 1:J) {",
                 "            for (k in 1:K) {", paste0(
                   "                log(phi[i, j, k]) <- ", term1, " + ", term2),
                 "            }",
                 "        }")
    }
  } else {
    if (phi == "i") {
      if (M_cov_phi == 1) {
        model <- c(model,
                   "        log(phi[i]) <- alpha[i, 1] * cov_phi[1]")
      } else {
        model <- c(model,
                   "        log(phi[i]) <- inprod(alpha[i, 1:M_phi], cov_phi[1:M_phi])")
      }
    } else if (phi == "ij") {
      if (M_cov_phi == 1) {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            log(phi[i, j]) <- alpha[i, 1] * cov_phi[j, 1]",
                   "        }")
      } else {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            log(phi[i, j]) <- inprod(alpha[i, 1:M_phi], cov_phi[j, 1:M_phi])",
                   "        }")
      }
    } else if (phi == "ijk") {
      if (M_cov_phi == 1) {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            for (k in 1:K) {",
                   "                log(phi[i, j, k]) <- alpha[i, 1] * cov_phi[j, k, 1]",
                   "            }",
                   "        }")
      } else {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            for (k in 1:K) {",
                   "                log(phi[i, j, k]) <- inprod(alpha[i, 1:M_phi], cov_phi[j, k, 1:M_phi])",
                   "            }",
                   "        }")
      }
    }
  }

  if (theta_shared) {
    if (theta == "i") {
      if (M_cov_theta == 1) {
        term1 <- "beta[i, 1] * cov_theta[1]"
      } else {
        term1 <- "inprod(beta[i, 1:M_theta], cov_theta[1:M_theta])"
      }
      if (M_cov_theta_shared == 1) {
        term2 <- "beta_shared[1] * cov_theta_shared[i, 1]"
      } else {
        term2 <- "inprod(beta_shared[1:M_theta_shared], cov_theta_shared[i, 1:M_theta_shared])"
      }
      model <- c(model, paste0(
        "        logit(theta[i]) <- ", term1, " + ", term2))
    } else if (theta == "ij") {
      if (M_cov_theta == 1) {
        term1 <- "beta[i, 1] * cov_theta[j, 1]"
      } else {
        term1 <- "inprod(beta[i, 1:M_theta], cov_theta[j, 1:M_theta])"
      }
      if (M_cov_theta_shared == 1) {
        term2 <- "beta_shared[1] * cov_theta_shared[i, j, 1]"
      } else {
        term2 <- "inprod(beta_shared[1:M_theta_shared], cov_theta_shared[i, j, 1:M_theta_shared])"
      }
      model <- c(model,
                 "        for (j in 1:J) {", paste0(
                   "            logit(theta[i, j]) <- ", term1, " + ", term2),
                 "        }")
    } else if (theta == "ijk") {
      if (M_cov_theta == 1) {
        term1 <- "beta[i, 1] * cov_theta[j, k, 1]"
      } else {
        term1 <- "inprod(beta[i, 1:M_theta], cov_theta[j, k, 1:M_theta])"
      }
      if (M_cov_theta_shared == 1) {
        term2 <- "beta_shared[1] * cov_theta_shared[i, j, k, 1]"
      } else {
        term2 <- "inprod(beta_shared[1:M_theta_shared], cov_theta_shared[i, j, k, 1:M_theta_shared])"
      }
      model <- c(model,
                 "        for (j in 1:J) {",
                 "            for (k in 1:K) {", paste0(
                   "                logit(theta[i, j, k]) <- ", term1, " + ", term2),
                 "            }",
                 "        }")
    }
  } else {
    if (theta == "i")
      if (M_cov_theta == 1) {
        model <- c(model,
                   "        logit(theta[i]) <- beta[i, 1] * cov_theta[1]")
      } else {
        model <- c(model,
                   "        logit(theta[i]) <- inprod(beta[i, 1:M_theta], cov_theta[1:M_theta])")
      }
    if (theta == "ij")
      if (M_cov_theta == 1) {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            logit(theta[i, j]) <- beta[i, 1] * cov_theta[j, 1]",
                   "        }")
      } else {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            logit(theta[i, j]) <- inprod(beta[i, 1:M_theta], cov_theta[j, 1:M_theta])",
                   "        }")
      }
    if (theta == "ijk")
      if (M_cov_theta == 1) {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            for (k in 1:K) {",
                   "                logit(theta[i, j, k]) <- beta[i, 1] * cov_theta[j, k, 1]",
                   "            }",
                   "        }")
      } else {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            for (k in 1:K) {",
                   "                logit(theta[i, j, k]) <- inprod(beta[i, 1:M_theta], cov_theta[j, k, 1:M_theta])",
                   "            }",
                   "        }")
      }
  }

  if (psi_shared) {
    if (psi == "i") {
      if (M_cov_psi == 1) {
        term1 <- "gamma[i, 1] * cov_psi[1]"
      } else {
        term1 <- "inprod(gamma[i, 1:M_psi], cov_psi[1:M_psi])"
      }
      if (M_cov_psi_shared == 1) {
        term2 <- "gamma_shared[1] * cov_psi_shared[i, 1]"
      } else {
        term2 <- "inprod(gamma_shared[1:M_psi_shared], cov_psi_shared[i, 1:M_psi_shared])"
      }
      model <- c(model,
                 paste0("        logit(psi[i]) <- ", term1, " + ", term2))

    } else if (psi == "ij") {
      if (M_cov_psi == 1) {
        term1 <- "gamma[i, 1] * cov_psi[j, 1]"
      } else {
        term1 <- "inprod(gamma[i, 1:M_psi], cov_psi[j, 1:M_psi])"
      }
      if (M_cov_psi_shared == 1) {
        term2 <- "gamma_shared[1] * cov_psi_shared[i, j, 1]"
      } else {
        term2 <- "inprod(gamma_shared[1:M_psi_shared], cov_psi_shared[i, j, 1:M_psi_shared])"
      }
      model <- c(model,
                 "        for (j in 1:J) {", paste0(
                   "            logit(psi[i, j]) <- ", term1, " + ", term2),
                 "        }")
    }
  } else {
    if (psi == "i") {
      if (M_cov_psi == 1) {
        model <- c(model,
                   "        logit(psi[i]) <- gamma[i, 1] * cov_psi[1]")
      } else {
        model <- c(model,
                   "        logit(psi[i]) <- inprod(gamma[i, 1:M_psi], cov_psi[1:M_psi])")
      }
    } else if (psi == "ij") {
      if (M_cov_psi == 1) {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            logit(psi[i, j]) <- gamma[i, 1] * cov_psi[j, 1]",
                   "        }")
      } else {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            logit(psi[i, j]) <- inprod(gamma[i, 1:M_psi], cov_psi[j, 1:M_psi])",
                   "        }")
      }
    }
  }

  model <- c(model,
             readLines(system.file("nimble",
                                   "occumb_template5.nimble",
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

attach_nimble_package <- function() {
  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop("Package 'nimble' is required. Please install it.", call. = FALSE)
  }

  is_attached <- ("package:nimble" %in% search())
  if (!is_attached) {
    attachNamespace("nimble")
  }
}

to_model_code <- function(model_code_strings) {
  model_code <- paste0(model_code_strings, collapse = "\n")
  model_code <- str2lang(model_code)
  model_code <- nimble::nimbleCode(model_code)
  model_code
}

set_const_nimble <- function(const, data) {
  # Original lengths before padding with a dummy value
  len_m_phi   <- length(data$m_phi)
  len_m_theta <- length(data$m_theta)
  len_m_psi   <- length(data$m_psi)

  const_nimble <- c(
    const[c("I", "J", "K", "N")],
    data[c("M", "M_phi_shared", "M_theta_shared", "M_psi_shared",
           "m_phi", "m_theta", "m_psi", "prior_prec", "prior_ulim")],
    M_phi = len_m_phi, M_theta = len_m_theta, M_psi = len_m_psi)
  const_nimble$m_phi   <- pad_dummy_value(data$m_phi)
  const_nimble$m_theta <- pad_dummy_value(data$m_theta)
  const_nimble$m_psi   <- pad_dummy_value(data$m_psi)
  const_nimble$rho_index <- make_rho_index(data$M)

  # Drop elements with NA names. These correspond to M_*_shared entries missing from data.
  const_nimble <- const_nimble[!is.na(names(const_nimble))]
  const_nimble
}

make_rho_index <- function(M) {
  index <- matrix(0L, M, M)
  ind <- 1L
  for (m1 in 1:(M - 1)) {
    for (m2 in (m1 + 1):M) {
      index[m1, m2] <- ind
      ind <- ind + 1L
    }
  }
  index
}

set_data_nimble <- function(data) {
  data_nimble <- data[c("y", "cov_phi", "cov_theta", "cov_psi",
                        "cov_phi_shared", "cov_theta_shared", "cov_psi_shared")]
  # Drop elements with NA names. These correspond to cov_*_shared entries missing from data.
  data_nimble <- data_nimble[!is.na(names(data_nimble))]
  data_nimble
}

pad_dummy_value <- function(x, dummy_value = -999) {
  if (length(x) == 1L) {
    c(x, dummy_value)
  } else {
    x
  }
}

set_inits_nimble <- function(inits, seed, n.chains, n_rho) {
  if (is.null(seed)) {
    seeds <- floor(stats::runif(n.chains, min = 1, max = 1e+05))
  } else if (isFALSE(seed) && length(seed) == 1L) {
    set.seed(NULL)
    seeds <- floor(stats::runif(n.chains, min = 1, max = 1e+05))
  } else if (isTRUE(seed) && length(seed) == 1L) {
    seeds <- seq_len(n.chains)
  } else if (is.numeric(seed) && length(seed) == n.chains) {
    if (length(unique(seed)) < n.chains) {
      nimble::messageIfVerbose("[Warning] 'seed' has duplicates; some chains may be identical.", call. = FALSE)
    }
    seeds <- seed
  } else if (is.numeric(seed) && length(seed) == 1L) {
    nimble::messageIfVerbose("[Note] Expanding a single 'seed' to per-chain seeds (seed, seed+1, ...).")
    seeds <- seed + seq_len(n.chains) - 1
  } else {
    stop("Invalid 'seed'. Use TRUE or a numeric vector of length n.chains. See nimble::runMCMC(setSeed = ...).", call. = FALSE)
  }
  inits_nimble <- lapply(seeds, function(s) {
    set.seed(s)
    i <- inits()
    i$rho <- double(n_rho)
    append(i, list(.RNG.seed = s))
  })
  inits_nimble
}


detect_cores_omit_one <- function() {
  # detect cores (may be NA on some systems)
  n <- parallel::detectCores()
  if (is.na(n)) stop("Cannot auto-detect n.cores; please set n.cores")
  # omit cores, but always keep at least one
  max(1L, n - 1L)
}

make_jagsui_compatible <- function(fit, env = parent.frame()) {
  with(env, {
    # Remove whitespace in parameter names (e.g., "alpha[1, 1]" -> "alpha[1,1]")
    rownames(fit$summary) <- gsub(pattern = "\\s", replacement = "",
                                  rownames(fit$summary))
    for (chain in seq_along(fit$samples)) {
      colnames(fit$samples[[chain]]) <- gsub(pattern = "\\s", replacement = "",
                                             colnames(fit$samples[[chain]]))
    }
    # Rename rho labels from vector style to matrix style
    # (e.g., "rho[1]" -> "rho[1,2]") to match JAGS convention
    rownames(fit$summary) <- rename_rho_labels(rownames(fit$summary),
                                               const_nimble$rho_index)
    for (chain in seq_along(fit$samples)) {
      colnames(fit$samples[[chain]]) <- rename_rho_labels(
        colnames(fit$samples[[chain]]), const_nimble$rho_index)
    }
    # Reshape sims.list$rho to [samples, M-1, M] to match JAGS convention
    fit$sims.list$rho <- reshape_rho_sims(fit$sims.list$rho,
                                          const_nimble$rho_index)
    # Strip "[1]" from scalar parameters to match JAGS convention
    # (e.g., "gamma_shared[1]" -> "gamma_shared")
    rownames(fit$summary) <- strip_scalar_indices(rownames(fit$summary))
    for (chain in seq_along(fit$samples)) {
      colnames(fit$samples[[chain]]) <-
        strip_scalar_indices(colnames(fit$samples[[chain]]))
    }
    for (nm in names(fit$sims.list)) {
      if (is.matrix(fit$sims.list[[nm]]) && ncol(fit$sims.list[[nm]]) == 1) {
        fit$sims.list[[nm]] <- drop(fit$sims.list[[nm]])
      }
    }
    fit$parallel   <- parallel
    fit$parameters <- params
    fit$model      <- to_occumb_nimble_model(model_code_strings, const_nimble, data_nimble)
    fit$modfile    <- model_file
    fit$run.date   <- start_time
    fit$mcmc.info$n.iter       <- n.iter
    fit$mcmc.info$n.burnin     <- n.burnin
    fit$mcmc.info$n.thin       <- n.thin
    fit$mcmc.info$elapsed.mins <- elapsed_mins
    if (store.data) {
      fit$data  <- data
      fit$inits <- inits_nimble
    }
    fit
  })
}

rename_rho_labels <- function(param_names, rho_index) {
  M <- nrow(rho_index)
  mapping <- character()
  for (m1 in 1:(M - 1)) {
    for (m2 in (m1 + 1):M) {
      k <- rho_index[m1, m2]
      mapping[paste0("rho[", k, "]")] <- paste0("rho[", m1, ",", m2, "]")
    }
  }
  idx <- param_names %in% names(mapping)
  param_names[idx] <- mapping[param_names[idx]]
  param_names
}

reshape_rho_sims <- function(rho_sims, rho_index) {
  M <- nrow(rho_index)
  n_draws <- nrow(rho_sims)
  out <- array(NA, dim = c(n_draws, M - 1, M))
  for (m1 in 1:(M - 1)) {
    for (m2 in (m1 + 1):M) {
      k <- rho_index[m1, m2]
      out[, m1, m2] <- rho_sims[, k]
    }
  }
  out
}

strip_scalar_indices <- function(param_names) {
  # Identify names like "foo[1]" where "foo[2]" does not exist
  scalar_pattern <- "^(.+)\\[1\\]$"
  candidate_indices <- grep(scalar_pattern, param_names)
  if (length(candidate_indices) == 0) return(param_names)

  base_names <- sub(scalar_pattern, "\\1", param_names[candidate_indices])
  is_scalar <- !(paste0(base_names, "[2]") %in% param_names)
  param_names[candidate_indices[is_scalar]] <- base_names[is_scalar]
  param_names
}

to_occumb_nimble_model <- function(model_code_strings, const, data) {
  structure(
    list(model_code_strings = model_code_strings, const = const, data = data),
    class = "occumb_nimble_model"
  )
}

#' @export
print.occumb_nimble_model <- function(x, ...) {
  cat(crayon::bold("NIMBLE model:"), "\n\n")
  for (i in seq_along(x$model_code_strings)) {
    cat(x$model_code_strings[i], "\n", sep = "")
  }
}
