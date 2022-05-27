additional_test <- FALSE

I <- 2 # Number of species
J <- 2 # Number of sites
K <- 2 # Number of replicates
data <- occumbData(
    y = array(sample.int(I * J * K), dim = c(I, J, K)),
    spec_cov = list(cov1 = rnorm(I)),
    site_cov = list(cov2 = rnorm(J),
                    cov3 = factor(1:J)),
    repl_cov = list(cov4 = matrix(rnorm(J * K), J, K)))

# Fitting a null model (includes only species-specific intercepts)
res0 <- occumb(data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)

# Add species-specific effects of site covariates in occupancy probabilities
res1 <- occumb(formula_psi = ~ cov2, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)
res1a <- occumb(formula_theta = ~ cov2, data = data,
                n.chains = 1, n.adapt = 0, n.burnin = 0,
                n.thin = 1, n.iter = 10, verbose = FALSE)
res1b <- occumb(formula_phi = ~ cov2, data = data,
                n.chains = 1, n.adapt = 0, n.burnin = 0,
                n.thin = 1, n.iter = 10, verbose = FALSE)
res2 <- occumb(formula_psi = ~ cov3, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)
res3 <- occumb(formula_psi = ~ cov2 * cov3, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)

# Add species covariate in the three parameters
res4 <- occumb(formula_phi_shared = ~ cov1, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)
res5 <- occumb(formula_theta_shared = ~ cov1, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)
res6 <- occumb(formula_psi_shared = ~ cov1, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)

# Add replicate covariates
res7 <- occumb(formula_phi = ~ cov4, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)
res8 <- occumb(formula_theta = ~ cov4, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)


### Tests for eval_util_L/R ----------------------------------------------------

test_that("eval_util_L() outputs a data frame with the additional Utility column", {
    settings <- data.frame(K = rep(1, 3), N = rep(1, 3), x = NA)

    # Null model
    test0 <- eval_util_L(settings, res0, cores = 1)
    checkmate::expect_data_frame(test0)
    expect_equal(colnames(test0), c(colnames(settings), "Utility"))
    expect_equal(test0[, -ncol(test0)], settings)

    # Model with site covariates
    test1 <- eval_util_L(settings, res1, cores = 1)
    checkmate::expect_data_frame(test1)
    expect_equal(colnames(test1), c(colnames(settings), "Utility"))
    expect_equal(test1[, -ncol(test1)], settings)

    test1a <- eval_util_L(settings, res1a, cores = 1)
    checkmate::expect_data_frame(test1a)
    expect_equal(colnames(test1a), c(colnames(settings), "Utility"))
    expect_equal(test1a[, -ncol(test1a)], settings)

    test1b <- eval_util_L(settings, res1b, cores = 1)
    checkmate::expect_data_frame(test1b)
    expect_equal(colnames(test1b), c(colnames(settings), "Utility"))
    expect_equal(test1b[, -ncol(test1b)], settings)

    test2 <- eval_util_L(settings, res2, cores = 1)
    checkmate::expect_data_frame(test2)
    expect_equal(colnames(test2), c(colnames(settings), "Utility"))
    expect_equal(test2[, -ncol(test2)], settings)

    test3 <- eval_util_L(settings, res3, cores = 1)
    checkmate::expect_data_frame(test3)
    expect_equal(colnames(test3), c(colnames(settings), "Utility"))
    expect_equal(test3[, -ncol(test3)], settings)

    # Model with species covariates
    test4 <- eval_util_L(settings, res4, cores = 1)
    checkmate::expect_data_frame(test4)
    expect_equal(colnames(test4), c(colnames(settings), "Utility"))
    expect_equal(test4[, -ncol(test4)], settings)

    test5 <- eval_util_L(settings, res5, cores = 1)
    checkmate::expect_data_frame(test5)
    expect_equal(colnames(test5), c(colnames(settings), "Utility"))
    expect_equal(test5[, -ncol(test5)], settings)

    test6 <- eval_util_L(settings, res6, cores = 1)
    checkmate::expect_data_frame(test6)
    expect_equal(colnames(test6), c(colnames(settings), "Utility"))
    expect_equal(test6[, -ncol(test6)], settings)
})

if (additional_test) {
    test_that("eval_util_R() outputs a data frame with the additional Utility column", {
        settings <- data.frame(J = rep(1, 3), K = rep(1, 3), N = rep(1, 3), x = NA)

        # Null model
        for (n in seq_len(1E4)) {
            test0 <- try(eval_util_R(settings, res0, cores = 1), silent = TRUE)
            if (class(test0) != "try-error") break
        }
        checkmate::expect_data_frame(test0)
        expect_equal(colnames(test0), c(colnames(settings), "Utility"))
        expect_equal(test0[, -ncol(test0)], settings)

        # Model with species covariates
        for (n in seq_len(1E4)) {
            test4 <- try(eval_util_R(settings, res4, cores = 1), silent = TRUE)
            if (class(test4) != "try-error") break
        }
        checkmate::expect_data_frame(test4)
        expect_equal(colnames(test4), c(colnames(settings), "Utility"))
        expect_equal(test4[, -ncol(test4)], settings)

        for (n in seq_len(1E4)) {
            test5 <- try(eval_util_R(settings, res5, cores = 1), silent = TRUE)
            if (class(test5) != "try-error") break
        }
        checkmate::expect_data_frame(test5)
        expect_equal(colnames(test5), c(colnames(settings), "Utility"))
        expect_equal(test5[, -ncol(test5)], settings)

        for (n in seq_len(1E4)) {
            test6 <- try(eval_util_R(settings, res6, cores = 1), silent = TRUE)
            if (class(test6) != "try-error") break
        }
        checkmate::expect_data_frame(test6)
        expect_equal(colnames(test6), c(colnames(settings), "Utility"))
        expect_equal(test6[, -ncol(test6)], settings)
    })
}

### Tests for qc_eval_util_L ---------------------------------------------------
test_that("qc_eval_util_L() blocks inappropriate settings", {
    expect_error(qc_eval_util_L(data.frame(Kx = rep(1, 2), N = rep(1, 2)), res0),
                 "The 'settings' argument does not contain column 'K'.")
    expect_error(qc_eval_util_L(data.frame(K = rep(1, 2), Nx = rep(1, 2)), res0),
                 "The 'settings' argument does not contain column 'N'.")
    expect_error(qc_eval_util_L(data.frame(K = rep(0, 2), N = rep(1, 2)), res0),
                 "'K' contains a non-positive value.")
    expect_error(qc_eval_util_L(data.frame(K = rep(1, 2), N = rep(0, 2)), res0),
                 "'N' contains a non-positive value.")
})

test_that("qc_eval_util_L() blocks models with replicate-specific parameters", {
    expect_error(qc_eval_util_L(data.frame(K = rep(1, 2), N = rep(1, 2)), res7),
                 "'phi' is replicate-specific: the current 'eval_util_L' is not applicable to models with replicate-specific parameters.")
    expect_error(qc_eval_util_L(data.frame(K = rep(1, 2), N = rep(1, 2)), res8),
                 "'theta' is replicate-specific: the current 'eval_util_L' is not applicable to models with replicate-specific parameters.")
})


### Tests for qc_eval_util_R ---------------------------------------------------
test_that("qc_eval_util_R() blocks inappropriate settings", {
    expect_error(qc_eval_util_R(data.frame(Jx = rep(1, 2),
                                           K = rep(1, 2),
                                           N = rep(1, 2)), res0),
                 "The 'settings' argument does not contain column 'J'.")
    expect_error(qc_eval_util_R(data.frame(J = rep(1, 2),
                                           Kx = rep(1, 2),
                                           N = rep(1, 2)), res0),
                 "The 'settings' argument does not contain column 'K'.")
    expect_error(qc_eval_util_R(data.frame(J = rep(1, 2),
                                           K = rep(1, 2),
                                           Nx = rep(1, 2)), res0),
                 "The 'settings' argument does not contain column 'N'.")
    expect_error(qc_eval_util_R(data.frame(J = rep(0, 2),
                                           K = rep(1, 2),
                                           N = rep(1, 2)), res0),
                 "'J' contains a non-positive value.")
    expect_error(qc_eval_util_R(data.frame(J = rep(1, 2),
                                           K = rep(0, 2),
                                           N = rep(1, 2)), res0),
                 "'K' contains a non-positive value.")
    expect_error(qc_eval_util_R(data.frame(J = rep(1, 2),
                                           K = rep(1, 2),
                                           N = rep(0, 2)), res0),
                 "'N' contains a non-positive value.")
})

test_that("qc_eval_util_R() blocks models with site-specific parameters", {
    expect_error(qc_eval_util_R(data.frame(J = rep(1, 2),
                                           K = rep(1, 2),
                                           N = rep(1, 2)), res1),
                 "'psi' is site-specific: the current 'eval_util_R' is not applicable to models with site-specific parameters.")
    expect_error(qc_eval_util_R(data.frame(J = rep(1, 2),
                                           K = rep(1, 2),
                                           N = rep(1, 2)), res1a),
                 "'theta' is site-specific: the current 'eval_util_R' is not applicable to models with site-specific parameters.")
    expect_error(qc_eval_util_R(data.frame(J = rep(1, 2),
                                           K = rep(1, 2),
                                           N = rep(1, 2)), res1b),
                 "'phi' is site-specific: the current 'eval_util_R' is not applicable to models with site-specific parameters.")
    expect_error(qc_eval_util_R(data.frame(J = rep(1, 2),
                                           K = rep(1, 2),
                                           N = rep(1, 2)), res2),
                 "'psi' is site-specific: the current 'eval_util_R' is not applicable to models with site-specific parameters.")
    expect_error(qc_eval_util_R(data.frame(J = rep(1, 2),
                                           K = rep(1, 2),
                                           N = rep(1, 2)), res3),
                 "'psi' is site-specific: the current 'eval_util_R' is not applicable to models with site-specific parameters.")
    expect_error(qc_eval_util_R(data.frame(J = rep(1, 2),
                                           K = rep(1, 2),
                                           N = rep(1, 2)), res7),
                 "'phi' is replicate-specific: the current 'eval_util_R' is not applicable to models with replicate-specific parameters.")
    expect_error(qc_eval_util_R(data.frame(J = rep(1, 2),
                                           K = rep(1, 2),
                                           N = rep(1, 2)), res8),
                 "'theta' is replicate-specific: the current 'eval_util_R' is not applicable to models with replicate-specific parameters.")
})


### Tests for eutil ------------------------------------------------------------
### * Tests are available only for non-parallel computations *

I <- 20; J <- 5; K <- 4; N <- 100; M <- 2
seed  <- rnorm(1)
z     <- array(rbinom(M * I * J, 1, 0.8), dim = c(M, I, J))
theta <- array(runif(M * I * J, min = 0.8), dim = c(M, I, J))
phi   <- array(rgamma(M * I * J, 1), dim = c(M, I, J))

test_that("eutil() works as expected for local scale", {
    # N_rep = 1
    ans <- with_seed(seed, {
        util_rep <- vector(length = M)
        for (n in seq_len(M))
            util_rep[n] <-
                cutil_local(z[n, , ], theta[n, , ], phi[n, , ], K, N)
        mean(util_rep)}
    )
    expect_equal(
        with_seed(seed, eutil(z, theta, phi, K, N, scale = "local",
                              N_rep = 1, cores = 1)),
        ans
    )
    # rep > 1
    ans <- with_seed(seed, {
        util_rep <- vector(length = M * 2)
        n <- rep(seq_len(M), each = 2)
        for (m in seq_along(n))
            util_rep[m] <-
                cutil_local(z[n[m], , ], theta[n[m], , ], phi[n[m], , ], K, N)
        mean(util_rep)}
    )
    expect_equal(
        with_seed(seed, eutil(z, theta, phi, K, N, scale = "local",
                              N_rep = 2, cores = 1)),
        ans
    )
})

test_that("eutil() works as expected for regional scale", {
    # rep = 1
    ans <- with_seed(seed, {
        util_rep <- vector(length = M)
        for (n in seq_len(M))
            util_rep[n] <-
                cutil_regional(z[n, , ], theta[n, , ], phi[n, , ], K, N)
        mean(util_rep)}
    )
    expect_equal(
        with_seed(seed, eutil(z, theta, phi, K, N, scale = "regional",
                              N_rep = 1, cores = 1)),
        ans
    )
    # rep > 1
    ans <- with_seed(seed, {
        util_rep <- vector(length = M * 2)
        n <- rep(seq_len(M), each = 2)
        for (m in seq_along(n))
            util_rep[m] <-
                cutil_regional(z[n[m], , ], theta[n[m], , ], phi[n[m], , ], K, N)
        mean(util_rep)}
    )
    expect_equal(
        with_seed(seed, eutil(z, theta, phi, K, N, scale = "regional",
                              N_rep = 2, cores = 1)),
        ans
    )
})

### Tests for cutil ------------------------------------------------------------

I <- 20; J <- 5; K <- 4; N <- 100
u <- r <- pi <- array(dim = c(I, J, K)); pi[1] <- NaN
while(any(is.nan(pi))) {
    z     <- matrix(rbinom(I * J, 1, 0.8), I, J)
    theta <- matrix(runif(I * J, min = 0.8), I, J)
    phi   <- matrix(rgamma(I * J, 1), I, J)
    seed  <- rnorm(1)
    set.seed(seed)

    for (k in seq_len(K)) {
        u[, , k] <- rbinom(I * J, 1, z * theta)
        r[, , k] <- rgamma(I * J, phi, 1)
    }
    for (k in seq_len(K)) {
        for (j in seq_len(J)) {
            for (i in seq_len(I)) {
                pi[i, j, k] <- u[i, j, k] * r[i, j, k] / sum(u[, j, k] * r[, j, k])
            }
        }
    }
}

test_that("cutil_local() works as expected", {
    ans <- with_seed(seed,
        sum(predict_detect_probs_local(predict_pi(z, theta, phi, K), N)) / J)
    expect_equal(with_seed(seed, cutil_local(z, theta, phi, K, N)), ans)
})

test_that("cutil_regional() works as expected", {
    ans <- with_seed(seed,
        sum(predict_detect_probs_regional(predict_pi(z, theta, phi, K), N)))
    expect_equal(with_seed(seed, cutil_regional(z, theta, phi, K, N)), ans)
})

test_that("predict_pi() works as expected", {
    with_seed(seed, expect_equal(predict_pi(z, theta, phi, K), pi))
    theta_zero <- matrix(0, I, J)
    expect_error(predict_pi(z, theta_zero, phi, K),
        "Failed to generate valid pi values under the given parameter set.")
})

test_that("predict_detect_probs_local() works as expected", {
    ans <- matrix(nrow = I, ncol = J)
    for (i in seq_len(I)) {
        for (j in seq_len(J)) {
            ans[i, j] <- 1 - prod((1 - pi[i, j, ])^N)
        }
    }
    expect_equal(predict_detect_probs_local(pi, N), ans)
})

test_that("predict_detect_probs_regional() works as expected", {
    ans <- vector(length = I)
    for (i in seq_len(I)) ans[i] <- 1 - prod((1 - pi[i, , ])^N)
    expect_equal(predict_detect_probs_regional(pi, N), ans)
})

