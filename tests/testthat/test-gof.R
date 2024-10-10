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
    out <- gof(fit, stats = "chi_squared", plot = FALSE, cores = 2)
    expect_identical(out@stats, "chi_squared")
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
    N_test  <- sum(y_test)
    expect_identical(Freeman_Tukey(y_test, N_test, pi_test),
                     sum(((sqrt(y_test)) - sqrt(N_test * pi_test))^2))
    expect_identical(chi_squared(y_test, N_test, pi_test),
                     sum((y_test - N_test * pi_test)^2 / (N_test * pi_test)))
})

test_that("Calculation of Bayes p-value is correct", {
    x <- rnorm(10); y <- rnorm(10)
    expect_identical(Bayesian_p_value(x, y), sum(x < y) / length(x))
})

### Tests for unbalanced designs ----------------------------------------------
data_unbalanced <- data
j_miss <- 1; k_miss <- 2
data_unbalanced@y[, j_miss, k_miss] <- 0
fit_unbalanced <- occumb(data = data_unbalanced,
                         n.chains = 1, n.burnin = 10, n.thin = 1, n.iter = 20,
                         verbose = FALSE)

test_that("Fit statistic is zero for missing observation", {
    y_ans  <- get_data(fit_unbalanced, "y")
    N_temp <- apply(y_ans, c(2, 3), sum)
    pi_temp <- get_post_samples(fit_unbalanced, "pi")

    test_Freeman_Tukey <- Freeman_Tukey(y_ans[, j_miss, k_miss],
                                        N_temp[j_miss, k_miss],
                                        pi_temp[1, , j_miss, k_miss])
    expect_identical(test_Freeman_Tukey, 0)
    test_deviance <- -2 * llmulti(y_ans[, j_miss, k_miss],
                                  N_temp[j_miss, k_miss],
                                  pi_temp[1, , j_miss, k_miss])
    expect_identical(test_deviance, 0)
    test_chi_squared <- chi_squared(y_ans[, j_miss, k_miss],
                                    N_temp[j_miss, k_miss],
                                    pi_temp[1, , j_miss, k_miss])
    expect_identical(test_chi_squared, 0)
})

test_that("get_y_rep() works with unbalanced data", {
    y_ans  <- get_data(fit_unbalanced, "y")
    I_temp <- dim(y_ans)[1]
    J_temp <- dim(y_ans)[2]
    K_temp <- dim(y_ans)[3]
    N_temp <- apply(y_ans, c(2, 3), sum)
    pi_temp <- get_post_samples(fit_unbalanced, "pi")
    M_temp  <- dim(pi_temp)[1]

    expect_no_error(get_y_rep(1, y_ans, N_temp, pi_temp))

    y_test <- get_y_rep(1, y_ans, N_temp, pi_temp)
    expect_identical(y_test[, j_miss, k_miss], y_ans[, j_miss, k_miss])
})

test_that("gof() works with unbalanced data", {
    expect_no_error(gof(fit_unbalanced, plot = FALSE))
})

