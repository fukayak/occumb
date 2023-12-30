### Test data ------------------------------------------------------------------
I <- 2; J <- 2; K <- 2
y <- array(sample.int(I * J * K), dim = c(I, J, K))
data  <- occumbData(y = y)
fit <- occumb(data = data,
              n.chains = 1, n.burnin = 10, n.thin = 1, n.iter = 20,
              verbose = FALSE)

### Tests for outputs ----------------------------------------------------------
test_that("Dimensions of the output are as expected", {
    out <- gof(fit, plot = FALSE, cores = 2)
    expect_identical(out@stats, "Freeman_Tukey")
    expect_equal(length(out@stats_obs), fit@fit$mcmc.info$n.samples)
    expect_equal(length(out@stats_rep), fit@fit$mcmc.info$n.samples)
    expect_false(any(c(out@p_value < 0, 1 < out@p_value)))
    out <- gof(fit, stats = "deviance", plot = FALSE, cores = 2)
    expect_identical(out@stats, "deviance")
    expect_equal(length(out@stats_obs), fit@fit$mcmc.info$n.samples)
    expect_equal(length(out@stats_rep), fit@fit$mcmc.info$n.samples)
    expect_false(any(c(out@p_value < 0, 1 < out@p_value)))
})

### Tests for quality controls -------------------------------------------------
test_that("Checks for data work", {
    expect_error(gof(fit = array(1, dim = rep(2, 3)), cores = 2),
                 "An occumbFit class object is expected for 'fit'")
})

### Tests for fit statistics ---------------------------------------------------
test_that("Calculation of fit statistics is correct", {
    y_test  <- sample.int(2)
    pi_test <- rep(0.5, 2)
    expect_identical(Freeman_Tukey(y_test, sum(y_test), pi_test),
                     sum(((sqrt(y_test)) - sqrt(sum(y_test) * pi_test))^2))
})

test_that("Calculation of Bayes p-value is correct", {
    x <- rnorm(10); y <- rnorm(10)
    expect_identical(Bayesian_p_value(x, y), sum(x < y) / length(x))
})

