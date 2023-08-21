### Tests for quality controls -------------------------------------------------
test_that("Checks for data work", {
    expect_error(occumb(data = array(1, dim = rep(2, 3))),
                 "An occumbData class object is expected for data")
})
test_that("Checks for formula work", {
    formulas <- c("formula_phi",
                  "formula_theta",
                  "formula_psi",
                  "formula_phi_shared",
                  "formula_theta_shared",
                  "formula_psi_shared")
    expect_error(occumb(1, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1,
                        data = occumbData(y = array(1, dim = rep(2, 3)))),
                 sprintf("Formula is expected for: %s", "formula_phi"))
    expect_error(occumb(~ 1, 1, ~ 1, ~ 1, ~ 1, ~ 1,
                        data = occumbData(y = array(1, dim = rep(2, 3)))),
                 sprintf("Formula is expected for: %s", "formula_theta"))
    expect_error(occumb(~ 1, ~ 1, 1, ~ 1, ~ 1, ~ 1,
                        data = occumbData(y = array(1, dim = rep(2, 3)))),
                 sprintf("Formula is expected for: %s", "formula_psi"))
    expect_error(occumb(~ 1, ~ 1, ~ 1, 1, ~ 1, ~ 1,
                        data = occumbData(y = array(1, dim = rep(2, 3)))),
                 sprintf("Formula is expected for: %s", "formula_phi_shared"))
    expect_error(occumb(~ 1, ~ 1, ~ 1, ~ 1, 1, ~ 1,
                        data = occumbData(y = array(1, dim = rep(2, 3)))),
                 sprintf("Formula is expected for: %s", "formula_theta_shared"))
    expect_error(occumb(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, 1,
                        data = occumbData(y = array(1, dim = rep(2, 3)))),
                 sprintf("Formula is expected for: %s", "formula_psi_shared"))
    expect_error(occumb(1, 1, 1, 1, 1, 1,
                        data = occumbData(y = array(1, dim = rep(2, 3)))),
                 sprintf("Formula is expected for: %s",
                         paste(formulas, collapse = ", ")))
})

### Tests for set_const() ------------------------------------------------------
test_that("Replacement of missing y to NA works", {
    y         <- array(1:8, dim = rep(2, 3))
    y[, 1, 1] <- 0
    result    <- set_const(occumbData(y = y))
    expect_identical(result$y[, 1, 1], as.numeric(rep(NA, 2)))
    expect_identical(result$y[, 2, 1], c(3, 4))
    expect_identical(result$y[, 1, 2], c(5, 6))
    expect_identical(result$y[, 2, 2], c(7, 8))
})
test_that("Replacement of sequence depth 0 to 1 works", {
    y         <- array(1:8, dim = rep(2, 3))
    y[, 1, 1] <- 0
    result    <- set_const(occumbData(y = y))
    expect_equal(result$N[1, 1], 1)
    expect_equal(result$N[2, 1], sum(3:4))
    expect_equal(result$N[1, 2], sum(5:6))
    expect_equal(result$N[2, 2], sum(7:8))
})

### Test cases for set_data/inits_code/set_params_monitored/write_jags_model ---
phi <- theta <- c("i", "ij", "ijk")
psi <- c("i", "ij")
phi_shared <- theta_shared <- psi_shared <- c(FALSE, TRUE)

cases <- expand.grid(phi, theta, psi, phi_shared, theta_shared, psi_shared)
colnames(cases) <- c("phi", "theta", "psi",
                     "phi_shared", "theta_shared", "psi_shared")

### Tests for set_data() -------------------------------------------------------
test_that("Data list is correct for 144 available models", {
    const <- list(I = sample.int(1E3, 1),
                  J = sample.int(1E3, 1),
                  K = sample.int(1E3, 1),
                  N = sample.int(1E3, 1),
                  y = sample.int(1E3, 1))
    prior_prec <- rnorm(1)
    prior_ulim <- rnorm(1)

    for (i in 1:nrow(cases)) {
        margs <- list(cov_phi          = rnorm(1),
                      cov_theta        = rnorm(1),
                      cov_psi          = rnorm(1),
                      cov_phi_shared   = rnorm(1),
                      cov_theta_shared = rnorm(1),
                      cov_psi_shared   = rnorm(1),
                      M                = sample.int(1E3, 1),
                      M_phi_shared     = sample.int(1E3, 1),
                      M_theta_shared   = sample.int(1E3, 1),
                      M_psi_shared     = sample.int(1E3, 1),
                      m_phi            = sample.int(1E3, 1),
                      m_theta          = sample.int(1E3, 1),
                      m_psi            = sample.int(1E3, 1),
                      phi_shared       = cases$phi_shared[i],
                      theta_shared     = cases$theta_shared[i],
                      psi_shared       = cases$psi_shared[i])

        ans <- list(I          = const$I,
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

        if (cases$phi_shared[i])
            ans$cov_phi_shared <- margs$cov_phi_shared
            ans$M_phi_shared   <- margs$M_phi_shared
        if (cases$theta_shared[i])
            ans$cov_theta_shared <- margs$cov_theta_shared
            ans$M_theta_shared   <- margs$M_theta_shared
        if (cases$psi_shared[i])
            ans$cov_psi_shared <- margs$cov_psi_shared
            ans$M_psi_shared   <- margs$M_psi_shared

        res <- set_data(const, margs, prior_prec, prior_ulim)
        expect_equal(res, ans)
    }
})

### Tests for inits_code() -----------------------------------------------------
test_that("Code for initial value function is correct for 144 available models", {
    const <- list(I = sample.int(1E3, 1),
                  J = sample.int(1E3, 1),
                  K = sample.int(1E3, 1),
                  N = sample.int(1E3, 1),
                  y = sample.int(1E3, 1))
    prior_prec <- rnorm(1)
    prior_ulim <- rnorm(1)

    for (i in 1:nrow(cases)) {
        margs <- list(cov_phi          = rnorm(1),
                      cov_theta        = rnorm(1),
                      cov_psi          = rnorm(1),
                      cov_phi_shared   = rnorm(1),
                      cov_theta_shared = rnorm(1),
                      cov_psi_shared   = rnorm(1),
                      M                = sample.int(1E3, 1),
                      M_phi_shared     = sample.int(1E3, 1),
                      M_theta_shared   = sample.int(1E3, 1),
                      M_psi_shared     = sample.int(1E3, 1),
                      m_phi            = sample.int(1E3, 1),
                      m_theta          = sample.int(1E3, 1),
                      m_psi            = sample.int(1E3, 1),
                      phi_shared       = cases$phi_shared[i],
                      theta_shared     = cases$theta_shared[i],
                      psi_shared       = cases$psi_shared[i])

        ans <- c(
            "function() {",
            "    list(z = matrix(1, const$I, const$J),",
            "         u = array(1, dim = c(const$I, const$J, const$K)),",
            "         r = array(stats::rnorm(const$I * const$J * const$K,",
            "                                mean = 1, sd = 0.1),",
            "                   dim = c(const$I, const$J, const$K)),")

        if (cases$phi_shared[i])
            ans <- c(ans, "         alpha_shared = stats::rnorm(margs$M_phi_shared, sd = 0.1),")
        if (cases$theta_shared[i])
            ans <- c(ans, "         beta_shared  = stats::rnorm(margs$M_theta_shared, sd = 0.1),")
        if (cases$psi_shared[i])
            ans <- c(ans, "         gamma_shared = stats::rnorm(margs$M_psi_shared, sd = 0.1),")

        ans <- c(ans,
            "         spec_eff = matrix(stats::rnorm(const$I * margs$M, sd = 0.1),",
            "                           const$I, margs$M),",
            "         Mu       = stats::rnorm(margs$M, sd = 0.1),",
            "         sigma    = stats::rnorm(margs$M, mean = 1, sd = 0.1))",
            "}")

        res <- inits_code(const, margs)
        expect_equal(res, ans)

    }
})

### Tests for set_params_monitored() -------------------------------------------
test_that("Parameter list is correct for 144 available models", {
    for (i in 1:nrow(cases)) {
        ans <- c("Mu", "sigma", "rho", "alpha", "beta", "gamma")

        if (cases$phi_shared[i])
            ans <- c(ans, "alpha_shared")
        if (cases$theta_shared[i])
            ans <- c(ans, "beta_shared")
        if (cases$psi_shared[i])
            ans <- c(ans, "gamma_shared")

        ans <- c(ans, c("phi", "theta", "psi", "z", "pi"))

        res <- set_params_monitored(phi_shared   = cases$phi_shared[i],
                                    theta_shared = cases$theta_shared[i],
                                    psi_shared   = cases$psi_shared[i])
        expect_equal(res, ans)
    }
})

### Tests for write_jags_model() -----------------------------------------------
test_that("JAGS code is correct for 144 available models", {
    for (i in 1:nrow(cases)) {
        ans <- readLines(system.file("jags",
                                     "occumb_template1.jags",
                                     package = "occumb"))

        if (cases$phi[i] == "i")
            ans <- c(ans,
                     "                r[i, j, k] ~ dgamma(phi[i], 1)")
        if (cases$phi[i] == "ij")
            ans <- c(ans,
                     "                r[i, j, k] ~ dgamma(phi[i, j], 1)")
        if (cases$phi[i] == "ijk")
            ans <- c(ans,
                     "                r[i, j, k] ~ dgamma(phi[i, j, k], 1)")

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template2.jags",
                                            package = "occumb")))

        if (cases$theta[i] == "i")
            ans <- c(ans,
                     "                u[i, j, k] ~ dbern(z[i, j] * theta[i])")
        if (cases$theta[i] == "ij")
            ans <- c(ans,
                     "                u[i, j, k] ~ dbern(z[i, j] * theta[i, j])")
        if (cases$theta[i] == "ijk")
            ans <- c(ans,
                     "                u[i, j, k] ~ dbern(z[i, j] * theta[i, j, k])")

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template3.jags",
                                            package = "occumb")))

        if (cases$psi[i] == "i")
            ans <- c(ans,
                     "            z[i, j] ~ dbern(psi[i])")
        if (cases$psi[i] == "ij")
            ans <- c(ans,
                     "            z[i, j] ~ dbern(psi[i, j])")

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template4.jags",
                                            package = "occumb")))

        if (cases$phi_shared[i]) {
            if (cases$phi[i] == "i")
                ans <- c(ans,
                         "        log(phi[i]) <- inprod(alpha[i, ], cov_phi[]) + inprod(alpha_shared[], cov_phi_shared[i, ])")
            if (cases$phi[i] == "ij")
                ans <- c(ans,
                         "        for (j in 1:J) {",
                         "            log(phi[i, j]) <- inprod(alpha[i, ], cov_phi[j, ]) + inprod(alpha_shared[], cov_phi_shared[i, j, ])",
                         "        }")
            if (cases$phi[i] == "ijk")
                ans <- c(ans,
                         "        for (j in 1:J) {",
                         "            for (k in 1:K) {",
                         "                log(phi[i, j, k]) <- inprod(alpha[i, ], cov_phi[j, k, ]) + inprod(alpha_shared[], cov_phi_shared[i, j, k, ])",
                         "            }",
                         "        }")
        } else {
            if (cases$phi[i] == "i")
                ans <- c(ans,
                         "        log(phi[i]) <- inprod(alpha[i, ], cov_phi[])")
            if (cases$phi[i] == "ij")
                ans <- c(ans,
                         "        for (j in 1:J) {",
                         "            log(phi[i, j]) <- inprod(alpha[i, ], cov_phi[j, ])",
                         "        }")
            if (cases$phi[i] == "ijk")
                ans <- c(ans,
                         "        for (j in 1:J) {",
                         "            for (k in 1:K) {",
                         "                log(phi[i, j, k]) <- inprod(alpha[i, ], cov_phi[j, k, ])",
                         "            }",
                         "        }")
        }

        if (cases$theta_shared[i]) {
            if (cases$theta[i] == "i")
                ans <- c(ans,
                         "        logit(theta[i]) <- inprod(beta[i, ], cov_theta[]) + inprod(beta_shared[], cov_theta_shared[i, ])")
            if (cases$theta[i] == "ij")
                ans <- c(ans,
                         "        for (j in 1:J) {",
                         "            logit(theta[i, j]) <- inprod(beta[i, ], cov_theta[j, ]) + inprod(beta_shared[], cov_theta_shared[i, j, ])",
                         "        }")
            if (cases$theta[i] == "ijk")
                ans <- c(ans,
                         "        for (j in 1:J) {",
                         "            for (k in 1:K) {",
                         "                logit(theta[i, j, k]) <- inprod(beta[i, ], cov_theta[j, k, ]) + inprod(beta_shared[], cov_theta_shared[i, j, k, ])",
                         "            }",
                         "        }")
        } else {
            if (cases$theta[i] == "i")
                ans <- c(ans,
                         "        logit(theta[i]) <- inprod(beta[i, ], cov_theta[])")
            if (cases$theta[i] == "ij")
                ans <- c(ans,
                         "        for (j in 1:J) {",
                         "            logit(theta[i, j]) <- inprod(beta[i, ], cov_theta[j, ])",
                         "        }")
            if (cases$theta[i] == "ijk")
                ans <- c(ans,
                         "        for (j in 1:J) {",
                         "            for (k in 1:K) {",
                         "                logit(theta[i, j, k]) <- inprod(beta[i, ], cov_theta[j, k, ])",
                         "            }",
                         "        }")
        }

        if (cases$psi_shared[i]) {
            if (cases$psi[i] == "i")
                ans <- c(ans,
                         "        logit(psi[i]) <- inprod(gamma[i, ], cov_psi[]) + inprod(gamma_shared[], cov_psi_shared[i, ])")
            if (cases$psi[i] == "ij")
                ans <- c(ans,
                         "        for (j in 1:J) {",
                         "            logit(psi[i, j]) <- inprod(gamma[i, ], cov_psi[j, ]) + inprod(gamma_shared[], cov_psi_shared[i, j, ])",
                         "        }")
        } else {
            if (cases$psi[i] == "i")
                ans <- c(ans,
                         "        logit(psi[i]) <- inprod(gamma[i, ], cov_psi[])")
            if (cases$psi[i] == "ij")
                ans <- c(ans,
                         "        for (j in 1:J) {",
                         "            logit(psi[i, j]) <- inprod(gamma[i, ], cov_psi[j, ])",
                         "        }")
        }

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template5.jags",
                                            package = "occumb")))

        if (cases$phi_shared[i])
            ans <- c(ans,
                     "    for (m in 1:M_phi_shared) {",
                     "        alpha_shared[m] ~ dnorm(0, prior_prec)",
                     "    }")
        if (cases$theta_shared[i])
            ans <- c(ans,
                     "    for (m in 1:M_theta_shared) {",
                     "        beta_shared[m] ~ dnorm(0, prior_prec)",
                     "    }")
        if (cases$psi_shared[i])
            ans <- c(ans,
                     "    for (m in 1:M_psi_shared) {",
                     "        gamma_shared[m] ~ dnorm(0, prior_prec)",
                     "    }")

        ans <- c(ans, "}", "")

        res <- write_jags_model(phi          = cases$phi[i],
                                theta        = cases$theta[i],
                                psi          = cases$psi[i],
                                phi_shared   = cases$phi_shared[i],
                                theta_shared = cases$theta_shared[i],
                                psi_shared   = cases$psi_shared[i])
        expect_equal(res, ans)
    }
})

### Test for get_data() --------------------------------------------------------
test_that("Outputs are correct when proper variable names are given", {
    I <- 2; J <- 2; K <- 2
    y <- array(sample.int(I * J * K), dim = c(I, J, K))
    spec_cov <- list(cov1 = rnorm(I))
    site_cov <- list(cov2 = rnorm(J), cov3 = factor(1:J))
    repl_cov <- list(cov4 = matrix(rnorm(J * K), J, K))
    data <- occumbData(
        y = y,
        spec_cov = spec_cov,
        site_cov = site_cov,
        repl_cov = repl_cov)

    fit <- occumb(data = data,
                  n.chains = 1, n.burnin = 10, n.thin = 1, n.iter = 20,
                  verbose = FALSE)

    expect_identical(get_data(fit, "y"), y)
    expect_identical(get_data(fit, "spec_cov"), spec_cov)
    expect_identical(get_data(fit, "site_cov"), site_cov)
    expect_identical(get_data(fit, "repl_cov"), repl_cov)
})

