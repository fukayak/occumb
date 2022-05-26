I <- 2 # Number of species
J <- 2 # Number of sites
K <- 2 # Number of replicates
data <- occumbData(
    y = array(sample.int(I * J * K), dim = c(I, J, K)),
    spec_cov = list(cov1 = rnorm(I)),
    site_cov = list(cov2 = rnorm(J),
                    cov3 = factor(1:J)),
    repl_cov = list(cov4 = matrix(rnorm(J * K), J, K)))

res0 <- occumb(data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)
res7 <- occumb(formula_phi = ~ cov4, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)
res8 <- occumb(formula_theta = ~ cov4, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)


### Tests for eval_util_L/R ----------------------------------------------------

test_that("eval_util_L() outputs a data frame with the additional Utility column", {
    settings <- data.frame(K = rep(1, 3), N = rep(1, 3), x = NA)
    test     <- eval_util_L(settings, res0)
    checkmate::expect_data_frame(test)
    expect_equal(colnames(test), c(colnames(settings), "Utility"))
    expect_equal(test[, -ncol(test)], settings)
})


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

