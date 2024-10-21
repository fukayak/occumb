set.seed(1)
### Test data ------------------------------------------------------------------
I <- 2
J <- 2
K <- 2
y <- array(sample.int(I * J * K), dim = c(I, J, K))
spec_cov <- list(cov1 = rnorm(I))
site_cov <- list(cov2 = rnorm(J), cov3 = factor(1:J))
repl_cov <- list(cov4 = matrix(rnorm(J * K), J, K))
data <- occumbData(
  y = y,
  spec_cov = spec_cov,
  site_cov = site_cov,
  repl_cov = repl_cov
)


### Test for occumbData summary() ----------------------------------------------
test_that("summary() works for occumbData as expected", {
  expect_snapshot(
    summary(data)
  )
})

### Test for occumbFit summary() -----------------------------------------------
test_that("summary() works for occumbFit as expected", {
  expect_snapshot(
    summary(occumb:::internal_fit)
  )
})
