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
              n.chains = 1, n.burnin = 10, n.thin = 1, n.iter = 20,
              verbose = FALSE)

### Tests for outputs ----------------------------------------------------------
test_that("Dimensions of the output are as expected", {
    out <- loglik(fit)
    expect_equal(nrow(out), fit@fit$mcmc.info$n.samples)
    expect_equal(ncol(out), dim(fit@data@y)[2] * dim(fit@data@y)[3])
})

### Tests for quality controls -------------------------------------------------
test_that("Checks for data work", {
    expect_error(loglik(fit = array(1, dim = rep(2, 3))),
                 "An occumbFit class object is expected for fit")
})

