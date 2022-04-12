### Tests for eutil ------------------------------------------------------------
I <- 4; J <- 3; K <- 2; N <- 100; M <- 2
seed  <- rnorm(1)
z     <- array(rbinom(M * I * J, 1, 0.8), dim = c(M, I, J))
theta <- array(runif(M * I * J, min = 0.5), dim = c(M, I, J))
phi   <- array(rgamma(M * I * J, 1), dim = c(M, I, J))

test_that("eutil() works as expected for local scale", {
    # rep = 1
    ans <- with_seed(seed, {
        util_rep <- vector(length = M)
        for (n in seq_len(M))
            util_rep[n] <-
                cutil_local(z[n, , ], theta[n, , ], phi[n, , ], K, N)
        mean(util_rep)})
    expect_equal(
        with_seed(seed, eutil(z, theta, phi, K, N, scale = "local")),
        ans)
    # rep > 1
    ans <- with_seed(seed, {
        util_rep <- vector(length = M * 2)
        n <- rep(seq_len(M), each = 2)
        for (m in seq_along(n))
            util_rep[m] <-
                cutil_local(z[n[m], , ], theta[n[m], , ], phi[n[m], , ], K, N)
        mean(util_rep)})
    expect_equal(
        with_seed(seed, eutil(z, theta, phi, K, N, scale = "local", rep = 2)),
        ans)
})

test_that("eutil() works as expected for regional scale", {
    # rep = 1
    ans <- with_seed(seed, {
        util_rep <- vector(length = M)
        for (n in seq_len(M))
            util_rep[n] <-
                cutil_regional(z[n, , ], theta[n, , ], phi[n, , ], K, N)
        mean(util_rep)})
    expect_equal(
        with_seed(seed, eutil(z, theta, phi, K, N, scale = "regional")),
        ans)
    # rep > 1
    ans <- with_seed(seed, {
        util_rep <- vector(length = M * 2)
        n <- rep(seq_len(M), each = 2)
        for (m in seq_along(n))
            util_rep[m] <-
                cutil_regional(z[n[m], , ], theta[n[m], , ], phi[n[m], , ], K, N)
        mean(util_rep)})
    expect_equal(
        with_seed(seed, eutil(z, theta, phi, K, N, scale = "regional", rep = 2)),
        ans)
})

### Tests for cutil ------------------------------------------------------------
I <- 4; J <- 3; K <- 2; N <- 100
u <- r <- pi <- array(dim = c(I, J, K)); pi[1] <- NaN
while(any(is.nan(pi))) {
    z     <- matrix(rbinom(I * J, 1, 0.8), I, J)
    theta <- matrix(runif(I * J, min = 0.5), I, J)
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

