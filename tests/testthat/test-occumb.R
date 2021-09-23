### Test for set_const() -------------------------------------------------------
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

### Test for set_modargs() -----------------------------------------------------
test_that("Setup for a null model works", {
    result <- set_modargs(~ 1, ~ 1, ~ 1,
                          occumbData(y = array(0, dim = rep(2, 3))))
    expect_equal(result$M, 3)
    expect_equal(result$code, "111000")
    expect_equal(result$cov_phi, 1)
    expect_equal(result$cov_theta, 1)
    expect_equal(result$cov_psi, 1)
    expect_equal(result$m_phi, 1)
    expect_equal(result$m_theta, 2)
    expect_equal(result$m_psi, 3)
})
test_that("Setup for varying-psi model works", {
    cov1 <- rnorm(2); cov2 <- factor(1:2)
    data <- occumbData(y = array(0, dim = rep(2, 3)),
                       site_cov = list(cov1 = cov1,
                                       cov2 = cov2))

    # Continuous effect
    result <- set_modargs(~ 1, ~ 1, ~ cov1, data)
    mm     <- model.matrix(as.formula(~ cov1), data@site_cov)
    expect_equal(result$M, 4)
    expect_equal(result$code, "112000")
    expect_equal(result$cov_phi, 1)
    expect_equal(result$cov_theta, 1)
    expect_equal(result$cov_psi, mm)
    expect_equal(result$m_phi, 1)
    expect_equal(result$m_theta, 2)
    expect_equal(result$m_psi, 3:4)

    # Discrete effect
    result <- set_modargs(~ 1, ~ 1, ~ cov2, data)
    mm     <- model.matrix(as.formula(~ cov2), data@site_cov)
    expect_equal(result$M, 4)
    expect_equal(result$code, "112000")
    expect_equal(result$cov_phi, 1)
    expect_equal(result$cov_theta, 1)
    expect_equal(result$cov_psi, mm)
    expect_equal(result$m_phi, 1)
    expect_equal(result$m_theta, 2)
    expect_equal(result$m_psi, 3:4)

    # Additive effects
    result <- set_modargs(~ 1, ~ 1, ~ cov1 + cov2, data)
    mm     <- model.matrix(as.formula(~ cov1 + cov2), data@site_cov)
    expect_equal(result$M, 5)
    expect_equal(result$code, "112000")
    expect_equal(result$cov_phi, 1)
    expect_equal(result$cov_theta, 1)
    expect_equal(result$cov_psi, mm)
    expect_equal(result$m_phi, 1)
    expect_equal(result$m_theta, 2)
    expect_equal(result$m_psi, 3:5)

    # Interactive effects
    result <- set_modargs(~ 1, ~ 1, ~ cov1 * cov2, data)
    mm     <- model.matrix(as.formula(~ cov1 * cov2), data@site_cov)
    expect_equal(result$M, 6)
    expect_equal(result$code, "112000")
    expect_equal(result$cov_phi, 1)
    expect_equal(result$cov_theta, 1)
    expect_equal(result$cov_psi, mm)
    expect_equal(result$m_phi, 1)
    expect_equal(result$m_theta, 2)
    expect_equal(result$m_psi, 3:6)
})
test_that("Setup for varying-theta/phi model works", {
    cov1 <- rnorm(2); cov2 <- factor(1:2)
    data <- occumbData(y = array(0, dim = rep(2, 3)),
                       site_cov = list(cov1 = cov1,
                                       cov2 = cov2))
    expect_error(set_modargs(~ 1, ~ cov1, ~ 1, data),
                 "Covariates in theta_formula are not yet supported")
    expect_error(set_modargs(~ cov1, ~ cov1, ~ 1, data),
                 "Covariates in phi_formula are not yet supported")
})

### Test for write_jags_model() ------------------------------------------------
phi   <- c("i", "ij", "ijk")
theta <- c("i", "ij", "ijk")
psi   <- c("i", "ij")

test_conditions <- expand.grid(phi, theta, psi)
colnames(test_conditions) <- c("phi", "theta", "psi")

test_that("JAGS codes are correct for models without shared effects", {
    for (i in 1:nrow(test_conditions)) {
        ans <- readLines(system.file("jags",
                                     "occumb_template1.jags",
                                     package = "occumb"))

        if (test_conditions$phi[i] == "i")
            ans <- c(ans,
                     "                x[i, j, k] ~ dgamma(phi[i], 1)")
        if (test_conditions$phi[i] == "ij")
            ans <- c(ans,
                     "                x[i, j, k] ~ dgamma(phi[i, j], 1)")
        if (test_conditions$phi[i] == "ijk")
            ans <- c(ans,
                     "                x[i, j, k] ~ dgamma(phi[i, j, k], 1)")

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template2.jags",
                                            package = "occumb")))

        if (test_conditions$theta[i] == "i")
            ans <- c(ans,
                     "                u[i, j, k] ~ dbern(z[i, j] * theta[i])")
        if (test_conditions$theta[i] == "ij")
            ans <- c(ans,
                     "                u[i, j, k] ~ dbern(z[i, j] * theta[i, j])")
        if (test_conditions$theta[i] == "ijk")
            ans <- c(ans,
                     "                u[i, j, k] ~ dbern(z[i, j] * theta[i, j, k])")

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template3.jags",
                                            package = "occumb")))

        if (test_conditions$psi[i] == "i")
            ans <- c(ans,
                     "            z[i, j] ~ dbern(psi[i])")
        if (test_conditions$psi[i] == "ij")
            ans <- c(ans,
                     "            z[i, j] ~ dbern(psi[i, j])")

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template4.jags",
                                            package = "occumb")))

        if (test_conditions$phi[i] == "i")
            ans <- c(ans,
                     "        log(phi[i]) <- inprod(alpha[i, ], cov_phi[])")
        if (test_conditions$phi[i] == "ij")
            ans <- c(ans,
                     "        log(phi[i, j]) <- inprod(alpha[i, ], cov_phi[i, j, ])")
        if (test_conditions$phi[i] == "ijk")
            ans <- c(ans,
                     "        log(phi[i, j, k]) <- inprod(alpha[i, ], cov_phi[i, j, k, ])")

        if (test_conditions$theta[i] == "i")
            ans <- c(ans,
                     "        logit(theta[i]) <- inprod(beta[i, ], cov_theta[])")
        if (test_conditions$theta[i] == "ij")
            ans <- c(ans,
                     "        logit(theta[i, j]) <- inprod(beta[i, ], cov_theta[i, j, ])")
        if (test_conditions$theta[i] == "ijk")
            ans <- c(ans,
                     "        logit(theta[i, j, k]) <- inprod(beta[i, ], cov_theta[i, j, k, ])")

        if (test_conditions$psi[i] == "i")
            ans <- c(ans,
                     "        logit(psi[i]) <- inprod(gamma[i, ], cov_psi[])")
        if (test_conditions$psi[i] == "ij")
            ans <- c(ans,
                     "        logit(psi[i, j]) <- inprod(gamma[i, ], cov_psi[i, j, ])")

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template5a.jags",
                                            package = "occumb")))

        res <- write_jags_model(phi = test_conditions$phi[i],
                                theta = test_conditions$theta[i],
                                psi = test_conditions$psi[i],
                                shared_effects = FALSE)
        expect_equal(res, ans)
    }
})
test_that("JAGS codes are correct for models with shared effects", {
    for (i in 1:nrow(test_conditions)) {
        ans <- readLines(system.file("jags",
                                     "occumb_template1.jags",
                                     package = "occumb"))

        if (test_conditions$phi[i] == "i")
            ans <- c(ans,
                     "                x[i, j, k] ~ dgamma(phi[i], 1)")
        if (test_conditions$phi[i] == "ij")
            ans <- c(ans,
                     "                x[i, j, k] ~ dgamma(phi[i, j], 1)")
        if (test_conditions$phi[i] == "ijk")
            ans <- c(ans,
                     "                x[i, j, k] ~ dgamma(phi[i, j, k], 1)")

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template2.jags",
                                            package = "occumb")))

        if (test_conditions$theta[i] == "i")
            ans <- c(ans,
                     "                u[i, j, k] ~ dbern(z[i, j] * theta[i])")
        if (test_conditions$theta[i] == "ij")
            ans <- c(ans,
                     "                u[i, j, k] ~ dbern(z[i, j] * theta[i, j])")
        if (test_conditions$theta[i] == "ijk")
            ans <- c(ans,
                     "                u[i, j, k] ~ dbern(z[i, j] * theta[i, j, k])")

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template3.jags",
                                            package = "occumb")))

        if (test_conditions$psi[i] == "i")
            ans <- c(ans,
                     "            z[i, j] ~ dbern(psi[i])")
        if (test_conditions$psi[i] == "ij")
            ans <- c(ans,
                     "            z[i, j] ~ dbern(psi[i, j])")

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template4.jags",
                                            package = "occumb")))

        if (test_conditions$phi[i] == "i")
            ans <- c(ans,
                     "        log(phi[i]) <- inprod(alpha[i, ], cov_phi[]) + inprod(alpha_shared[], cov_phi_shared[i, ])")
        if (test_conditions$phi[i] == "ij")
            ans <- c(ans,
                     "        log(phi[i, j]) <- inprod(alpha[i, ], cov_phi[i, j, ]) + inprod(alpha_shared[], cov_phi_shared[i, j, ])")
        if (test_conditions$phi[i] == "ijk")
            ans <- c(ans,
                     "        log(phi[i, j, k]) <- inprod(alpha[i, ], cov_phi[i, j, k, ]) + inprod(alpha_shared[], cov_phi_shared[i, j, k, ])")

        if (test_conditions$theta[i] == "i")
            ans <- c(ans,
                     "        logit(theta[i]) <- inprod(beta[i, ], cov_theta[]) + inprod(beta_shared[], cov_theta_shared[i, ])")
        if (test_conditions$theta[i] == "ij")
            ans <- c(ans,
                     "        logit(theta[i, j]) <- inprod(beta[i, ], cov_theta[i, j, ]) + inprod(beta_shared[], cov_theta_shared[i, j, ])")
        if (test_conditions$theta[i] == "ijk")
            ans <- c(ans,
                     "        logit(theta[i, j, k]) <- inprod(beta[i, ], cov_theta[i, j, k, ]) + inprod(beta_shared[], cov_theta_shared[i, j, k, ])")

        if (test_conditions$psi[i] == "i")
            ans <- c(ans,
                     "        logit(psi[i]) <- inprod(gamma[i, ], cov_psi[]) + inprod(gamma_shared[], cov_psi_shared[i, ])")
        if (test_conditions$psi[i] == "ij")
            ans <- c(ans,
                     "        logit(psi[i, j]) <- inprod(gamma[i, ], cov_psi[i, j, ]) + inprod(gamma_shared[], cov_psi_shared[i, j, ])")

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template5b.jags",
                                            package = "occumb")))

        res <- write_jags_model(phi = test_conditions$phi[i],
                                theta = test_conditions$theta[i],
                                psi = test_conditions$psi[i],
                                shared_effects = TRUE)
        expect_equal(res, ans)
    }
})

