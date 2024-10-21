## code to prepare `snapshot_occumb` object goes here
setwd("<occumb directory>")
devtools::load_all()

set.seed(1)
### Test data
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

internal_fit <- occumb(
  data = data,
  n.chains = 1, n.burnin = 10, n.thin = 1, n.iter = 20,
  verbose = FALSE
)

usethis::use_data(internal_fit, internal = TRUE)
