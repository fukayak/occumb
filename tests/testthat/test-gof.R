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
    out <- gof(fit, data, plot = FALSE)
    expect_identical(names(out), c("p_values", "stats_obs", "stats_rep"))
    expect_identical(names(out$p_values), c("deviance", "Freeman_Tukey"))
    expect_identical(names(out$stats_obs), c("deviance", "Freeman_Tukey"))
    expect_identical(names(out$stats_rep), c("deviance", "Freeman_Tukey"))
    expect_equal(length(out$stats_obs$deviance), fit@fit$mcmc.info$n.samples)
    expect_equal(length(out$stats_obs$Freeman_Tukey), fit@fit$mcmc.info$n.samples)
    expect_equal(length(out$stats_rep$deviance), fit@fit$mcmc.info$n.samples)
    expect_equal(length(out$stats_rep$Freeman_Tukey), fit@fit$mcmc.info$n.samples)
    expect_false(any(out$p_values < 0 & 1 < out$p_values))
})

### Tests for quality controls -------------------------------------------------
test_that("Checks for data work", {
    expect_error(gof(fit = array(1, dim = rep(2, 3)), data = data),
                 "An occumbFit class object is expected for fit")
    expect_error(gof(fit, data = array(1, dim = rep(2, 3))),
                 "An occumbData class object is expected for data")
    expect_error(gof(fit, data = data2),
                 "Dimension mismatch between the data and posterior: make sure to supply the dataset to which the model was applied")
})

### Tests for fit statistics ---------------------------------------------------
test_that("Calculations of fit statistics are correct", {
    y_test  <- sample.int(2)
    pi_test <- rep(0.5, 2)
    expect_identical(Freeman_Tukey(y_test, sum(y_test), pi_test),
                     sum(((sqrt(y_test)) - sqrt(sum(y_test) * pi_test))^2))
    expect_identical(llmulti(y_test, sum(y_test), pi_test),
                     dmultinom(y_test, sum(y_test), pi_test, log = TRUE))
})

test_that("Calculation of Bayes p-value is correct", {
    x <- rnorm(10); y <- rnorm(10)
    expect_identical(Bayesian_p_value(x, y), sum(x < y) / length(x))
})

