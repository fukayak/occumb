### Test data ------------------------------------------------------------------
I <- 2; J <- 2; K <- 2
y <- array(sample.int(I * J * K), dim = c(I, J, K))
spec_cov <- list(cov1 = rnorm(I))
site_cov <- list(cov2 = rnorm(J), cov3 = factor(1:J))
repl_cov <- list(cov4 = matrix(rnorm(J * K), J, K))
data <- occumbData(y = y,
                   spec_cov = spec_cov,
                   site_cov = site_cov,
                   repl_cov = repl_cov)
data2 <- occumbData(y = array(sample.int((I + 1) * J * K), dim = c(I + 1, J, K)),
                    spec_cov = list(cov1 = rnorm(I + 1)),
                    site_cov = site_cov,
                    repl_cov = repl_cov)

fit <- occumb(data = data,
              n.chains = 2, n.burnin = 10, n.thin = 1, n.iter = 20,
              verbose = FALSE)

### Tests for the MCMC output sorting ------------------------------------------
test_that("The pi array is correctly aligned", {
    pi <- get_post_samples(fit, "pi")

    for (i in seq_len(I)) {
        for (j in seq_len(J)) {
            for (k in seq_len(K)) {
                ans <- c(fit@fit$samples[[1]][, sprintf("pi[%s,%s,%s]", i, j, k)],
                         fit@fit$samples[[2]][, sprintf("pi[%s,%s,%s]", i, j, k)])
                # Order: samples in chain1, samples in chain2, ... 
                expect_identical(pi[, i, j, k], ans)
            }
        }
    }
})

test_that("get_m() works correctly", {
    pi <- get_post_samples(fit, "pi")
    n_iter  <- fit@fit$mcmc.info$n.samples / fit@fit$mcmc.info$n.chains
    n_chain <- fit@fit$mcmc.info$n.chains

    for (i in seq_len(I)) {
        for (j in seq_len(J)) {
            for (k in seq_len(K)) {
                ans <- c(fit@fit$samples[[1]][, sprintf("pi[%s,%s,%s]", i, j, k)],
                         fit@fit$samples[[2]][, sprintf("pi[%s,%s,%s]", i, j, k)])
                for (iter in seq_len(n_iter)) {
                    for (chain in seq_len(n_chain)) {
                        m <- get_m(iter, chain, n_iter)
                        expect_identical(pi[m, i, j, k], ans[m])
                    }
                }
            }
        }
    }
})

### Tests for outputs ----------------------------------------------------------
test_that("Dimensions of the output are as expected", {
    out <- loglik(fit)
    expect_equal(dim(out)[1], fit@fit$mcmc.info$n.samples / fit@fit$mcmc.info$n.chains)
    expect_equal(dim(out)[2], fit@fit$mcmc.info$n.chains)
    expect_equal(dim(out)[3], dim(fit@data@y)[2] * dim(fit@data@y)[3])
})

### Tests for quality controls -------------------------------------------------
test_that("Checks for data work", {
    expect_error(loglik(fit = array(1, dim = rep(2, 3))),
                 "An occumbFit class object is expected for fit")
})

