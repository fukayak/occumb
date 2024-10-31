### Test data ------------------------------------------------------------------
I <- 2
J <- 2
K <- 2
y <- array(sample.int(I * J * K), dim = c(I, J, K))
spec_cov <- list(cov1 = rnorm(I))
site_cov <- list(cov2 = rnorm(J), cov3 = factor(1:J))
repl_cov <- list(cov4 = matrix(rnorm(J * K), J, K))
data <- occumbData(y = y,
                   spec_cov = spec_cov,
                   site_cov = site_cov,
                   repl_cov = repl_cov)

fit <- occumb(data = data,
              n.chains = 1, n.burnin = 10, n.thin = 1, n.iter = 20,
              verbose = FALSE)

### Tests for assert_occumbFit() -----------------------------------------------
test_that("assert_occumbFit() works correctly", {
  expect_invisible(assert_occumbFit(fit))
  expect_error(assert_occumbFit(fit = array(1, dim = rep(2, 3))),
               "An occumbFit class object is expected for 'fit'")
})

### Tests for llmulti() --------------------------------------------------------
test_that("Calculation of multinomial log-likelihood is correct", {
  y_test  <- sample.int(2)
  pi_test <- rep(0.5, 2)
  expect_identical(llmulti(y_test, sum(y_test), pi_test),
                   dmultinom(y_test, sum(y_test), pi_test, log = TRUE))
})
