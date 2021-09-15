### Test for set_const() --------------------------------------------
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

### Test for set_modargs() --------------------------------------------
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

