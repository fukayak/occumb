### Test data ------------------------------------------------------------------
I <- 2
J <- 3
K <- 4
y <- array(sample.int(I * J * K), dim = c(I, J, K))
dimnames(y)[[1]] <- sprintf("species %s", 1:I)
data <- occumbData(y = y,
                   spec_cov = list(cov1 = rnorm(I),
                                   cov2 = letters[1:I],
                                   cov3 = factor(1:I),
                                   cov4 = c(TRUE, FALSE)),
                   site_cov = list(cov5 = rnorm(J),
                                   cov6 = letters[1:J],
                                   cov7 = factor(1:J),
                                   cov8 = c(TRUE, FALSE, TRUE)),
                   repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                                   cov10 = matrix(letters[1:(J * K)], J, K),
                                   cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))

### Tests for check_newdata() --------------------------------------------------
test_that("check_newdata() works correctly", {
  fit <- occumb(data = data,
                n.chains = 1, n.burnin = 10, n.thin = 1, n.iter = 20,
                verbose = FALSE)

  expect_no_error(check_newdata(data, fit))

  ## Check for mismatches in the number of species
  newdata <-
    occumbData(y = array(sample.int((I + 1) * J * K), dim = c(I + 1, J, K)),
               spec_cov = list(cov1 = rnorm(I + 1),
                               cov2 = letters[1:(I + 1)],
                               cov3 = factor(1:(I + 1)),
                               cov4 = c(TRUE, FALSE, FALSE)),
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J],
                               cov7 = factor(1:J),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit),
               "The number of species in 'newdata' \\(3\\) differs from that in the fitted data \\(2\\).")

  ## Check for covariate names and their order
  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1x = rnorm(I),
                               cov2x = letters[1:I],
                               cov3x = factor(1:I),
                               cov4x = c(TRUE, FALSE)),
               site_cov = list(cov5x = rnorm(J),
                               cov6x = letters[1:J],
                               cov7x = factor(1:J),
                               cov8x = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9x = matrix(rnorm(J * K), J, K),
                               cov10x = matrix(letters[1:(J * K)], J, K),
                               cov11x = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov1x, cov2x, cov3x, cov4x \\(newdata\\); cov1, cov2, cov3, cov4 \\(fitted data\\)
  site_cov: cov5x, cov6x, cov7x, cov8x \\(newdata\\); cov5, cov6, cov7, cov8 \\(fitted data\\)
  repl_cov: cov9x, cov10x, cov11x \\(newdata\\); cov9, cov10, cov11 \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov2 = rnorm(I),
                               cov1 = letters[1:I],
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov6 = rnorm(J),
                               cov5 = letters[1:J],
                               cov7 = factor(1:J),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov2, cov1, cov3, cov4 \\(newdata\\); cov1, cov2, cov3, cov4 \\(fitted data\\)
  site_cov: cov6, cov5, cov7, cov8 \\(newdata\\); cov5, cov6, cov7, cov8 \\(fitted data\\)")

  ## Check for covariate classes
  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = letters[1:I],
                               cov2 = rnorm(I),
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J],
                               cov7 = c(TRUE, FALSE, TRUE),
                               cov8 = factor(1:J)),
               repl_cov = list(cov9 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K),
                               cov10 = matrix(rnorm(J * K), J, K),
                               cov11 = matrix(letters[1:(J * K)], J, K)))
  expect_error(check_newdata(newdata, fit),
               "The covariate classes in 'newdata' must match those in the fitted data.
  cov1: character \\(newdata\\), numeric \\(fitted data\\)
  cov2: numeric \\(newdata\\), character \\(fitted data\\)
  cov7: logical \\(newdata\\), factor \\(fitted data\\)
  cov8: factor \\(newdata\\), logical \\(fitted data\\)
  cov9: logical \\(newdata\\), numeric \\(fitted data\\)
  cov10: numeric \\(newdata\\), character \\(fitted data\\)
  cov11: character \\(newdata\\), logical \\(fitted data\\)")

  ## Check for new levels in discrete covariates
  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = rnorm(I),
                               cov2 = letters[1:I + 1],
                               cov3 = factor(1:I + 1),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J + 1],
                               cov7 = factor(1:J + 1),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K) + 1], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit),
               "The levels of discrete covariates in 'newdata' must match those in the fitted data.
  cov2: b, c \\(newdata\\); a, b \\(fitted data\\)
  cov3: 2, 3 \\(newdata\\); 1, 2 \\(fitted data\\)
  cov6: b, c, d \\(newdata\\); a, b, c \\(fitted data\\)
  cov7: 2, 3, 4 \\(newdata\\); 1, 2, 3 \\(fitted data\\)
  cov10: b, c, d, e, f, g, h, i, j, k, l, m \\(newdata\\); a, b, c, d, e, f, g, h, i, j, k, l \\(fitted data\\)")

  ## Check for orders of factor levels
  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = rnorm(I),
                               cov2 = letters[1:I],
                               cov3 = factor(1:I, levels = I:1),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J],
                               cov7 = factor(1:J, levels = J:1),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit),
               "The levels of discrete covariates in 'newdata' must match those in the fitted data.
  cov3: 1, 2 \\(newdata\\); 1, 2 \\(fitted data\\)
  cov7: 1, 2, 3 \\(newdata\\); 1, 2, 3 \\(fitted data\\)")

  ## Check for mismatches in the list of species names
  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = rnorm(I),
                               cov2 = letters[1:I],
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J],
                               cov7 = factor(1:J),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_warning(check_newdata(newdata, fit),
                 "The list of species names in 'newdata' does not match that in the fitted data; the list of species names in the fitted data will be added to the 'label' attribute of the returned object.")


  ### Additional tests for check for covariate names and their order
  ## Without spec_cov
  data_no_spec_cov <-
    occumbData(y = y,
               spec_cov = NULL,
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J],
                               cov7 = factor(1:J),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))

  fit_no_spec_cov <-
    occumb(data = data_no_spec_cov,
           n.chains = 1, n.burnin = 10, n.thin = 1, n.iter = 20,
           verbose = FALSE)

  expect_no_error(check_newdata(data_no_spec_cov, fit_no_spec_cov))

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = NULL,
               site_cov = list(cov5x = rnorm(J),
                               cov6x = letters[1:J],
                               cov7x = factor(1:J),
                               cov8x = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9x = matrix(rnorm(J * K), J, K),
                               cov10x = matrix(letters[1:(J * K)], J, K),
                               cov11x = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit_no_spec_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  site_cov: cov5x, cov6x, cov7x, cov8x \\(newdata\\); cov5, cov6, cov7, cov8 \\(fitted data\\)
  repl_cov: cov9x, cov10x, cov11x \\(newdata\\); cov9, cov10, cov11 \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = NULL,
               site_cov = list(cov6 = rnorm(J),
                               cov5 = letters[1:J],
                               cov7 = factor(1:J),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit_no_spec_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  site_cov: cov6, cov5, cov7, cov8 \\(newdata\\); cov5, cov6, cov7, cov8 \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = rnorm(I)),
               site_cov = list(cov5x = rnorm(J),
                               cov6x = letters[1:J],
                               cov7x = factor(1:J),
                               cov8x = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9x = matrix(rnorm(J * K), J, K),
                               cov10x = matrix(letters[1:(J * K)], J, K),
                               cov11x = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit_no_spec_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov1 \\(newdata\\); \\(None\\) \\(fitted data\\)
  site_cov: cov5x, cov6x, cov7x, cov8x \\(newdata\\); cov5, cov6, cov7, cov8 \\(fitted data\\)
  repl_cov: cov9x, cov10x, cov11x \\(newdata\\); cov9, cov10, cov11 \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = rnorm(I)),
               site_cov = list(cov6 = rnorm(J),
                               cov5 = letters[1:J],
                               cov7 = factor(1:J),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit_no_spec_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov1 \\(newdata\\); \\(None\\) \\(fitted data\\)
  site_cov: cov6, cov5, cov7, cov8 \\(newdata\\); cov5, cov6, cov7, cov8 \\(fitted data\\)")

  # Check for covariate classes
  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = NULL,
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J],
                               cov7 = c(TRUE, FALSE, TRUE),
                               cov8 = factor(1:J)),
               repl_cov = list(cov9 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K),
                               cov10 = matrix(rnorm(J * K), J, K),
                               cov11 = matrix(letters[1:(J * K)], J, K)))
  expect_error(check_newdata(newdata, fit_no_spec_cov),
               "The covariate classes in 'newdata' must match those in the fitted data.
  cov7: logical \\(newdata\\), factor \\(fitted data\\)
  cov8: factor \\(newdata\\), logical \\(fitted data\\)
  cov9: logical \\(newdata\\), numeric \\(fitted data\\)
  cov10: numeric \\(newdata\\), character \\(fitted data\\)
  cov11: character \\(newdata\\), logical \\(fitted data\\)")

  # Check for new levels in discrete covariates
  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = NULL,
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J + 1],
                               cov7 = factor(1:J + 1),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K) + 1], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit_no_spec_cov),
               "The levels of discrete covariates in 'newdata' must match those in the fitted data.
  cov6: b, c, d \\(newdata\\); a, b, c \\(fitted data\\)
  cov7: 2, 3, 4 \\(newdata\\); 1, 2, 3 \\(fitted data\\)
  cov10: b, c, d, e, f, g, h, i, j, k, l, m \\(newdata\\); a, b, c, d, e, f, g, h, i, j, k, l \\(fitted data\\)")

  # Check for orders of factor levels
  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = NULL,
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J],
                               cov7 = factor(1:J, levels = J:1),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit_no_spec_cov),
               "The levels of discrete covariates in 'newdata' must match those in the fitted data.
  cov7: 1, 2, 3 \\(newdata\\); 1, 2, 3 \\(fitted data\\)")

  ## Without site_cov
  data_no_site_cov <-
    occumbData(y = y,
               spec_cov = list(cov1 = rnorm(I),
                               cov2 = letters[1:I],
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = NULL,
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))

  fit_no_site_cov <-
    occumb(data = data_no_site_cov,
           n.chains = 1, n.burnin = 10, n.thin = 1, n.iter = 20,
           verbose = FALSE)

  expect_no_error(check_newdata(data_no_site_cov, fit_no_site_cov))

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1x = rnorm(I),
                               cov2x = letters[1:I],
                               cov3x = factor(1:I),
                               cov4x = c(TRUE, FALSE)),
               site_cov = NULL,
               repl_cov = list(cov9x = matrix(rnorm(J * K), J, K),
                               cov10x = matrix(letters[1:(J * K)], J, K),
                               cov11x = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit_no_site_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov1x, cov2x, cov3x, cov4x \\(newdata\\); cov1, cov2, cov3, cov4 \\(fitted data\\)
  repl_cov: cov9x, cov10x, cov11x \\(newdata\\); cov9, cov10, cov11 \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov2 = rnorm(I),
                               cov1 = letters[1:I],
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = NULL,
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit_no_site_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov2, cov1, cov3, cov4 \\(newdata\\); cov1, cov2, cov3, cov4 \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1x = rnorm(I),
                               cov2x = letters[1:I],
                               cov3x = factor(1:I),
                               cov4x = c(TRUE, FALSE)),
               site_cov = list(cov5 = rnorm(J)),
               repl_cov = list(cov9x = matrix(rnorm(J * K), J, K),
                               cov10x = matrix(letters[1:(J * K)], J, K),
                               cov11x = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit_no_site_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov1x, cov2x, cov3x, cov4x \\(newdata\\); cov1, cov2, cov3, cov4 \\(fitted data\\)
  site_cov: cov5 \\(newdata\\); \\(None\\) \\(fitted data\\)
  repl_cov: cov9x, cov10x, cov11x \\(newdata\\); cov9, cov10, cov11 \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov2 = rnorm(I),
                               cov1 = letters[1:I],
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov5 = rnorm(J)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit_no_site_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov2, cov1, cov3, cov4 \\(newdata\\); cov1, cov2, cov3, cov4 \\(fitted data\\)
  site_cov: cov5 \\(newdata\\); \\(None\\) \\(fitted data\\)")

  # Check for covariate classes
  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = letters[1:I],
                               cov2 = rnorm(I),
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = NULL,
               repl_cov = list(cov9 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K),
                               cov10 = matrix(rnorm(J * K), J, K),
                               cov11 = matrix(letters[1:(J * K)], J, K)))
  expect_error(check_newdata(newdata, fit_no_site_cov),
               "The covariate classes in 'newdata' must match those in the fitted data.
  cov1: character \\(newdata\\), numeric \\(fitted data\\)
  cov2: numeric \\(newdata\\), character \\(fitted data\\)
  cov9: logical \\(newdata\\), numeric \\(fitted data\\)
  cov10: numeric \\(newdata\\), character \\(fitted data\\)
  cov11: character \\(newdata\\), logical \\(fitted data\\)")

  # Check for new levels in discrete covariates
  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = rnorm(I),
                               cov2 = letters[1:I + 1],
                               cov3 = factor(1:I + 1),
                               cov4 = c(TRUE, FALSE)),
               site_cov = NULL,
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K) + 1], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit_no_site_cov),
               "The levels of discrete covariates in 'newdata' must match those in the fitted data.
  cov2: b, c \\(newdata\\); a, b \\(fitted data\\)
  cov3: 2, 3 \\(newdata\\); 1, 2 \\(fitted data\\)
  cov10: b, c, d, e, f, g, h, i, j, k, l, m \\(newdata\\); a, b, c, d, e, f, g, h, i, j, k, l \\(fitted data\\)")

  # Check for orders of factor levels
  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = rnorm(I),
                               cov2 = letters[1:I],
                               cov3 = factor(1:I, levels = I:1),
                               cov4 = c(TRUE, FALSE)),
               site_cov = NULL,
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit_no_site_cov),
               "The levels of discrete covariates in 'newdata' must match those in the fitted data.
  cov3: 1, 2 \\(newdata\\); 1, 2 \\(fitted data\\)")

  ## Without repl_cov
  data_no_repl_cov <-
    occumbData(y = y,
               spec_cov = list(cov1 = rnorm(I),
                               cov2 = letters[1:I],
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J],
                               cov7 = factor(1:J),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = NULL)

  fit_no_repl_cov <-
    occumb(data = data_no_repl_cov,
           n.chains = 1, n.burnin = 10, n.thin = 1, n.iter = 20,
           verbose = FALSE)

  expect_no_error(check_newdata(data_no_repl_cov, fit_no_repl_cov))

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1x = rnorm(I),
                               cov2x = letters[1:I],
                               cov3x = factor(1:I),
                               cov4x = c(TRUE, FALSE)),
               site_cov = list(cov5x = rnorm(J),
                               cov6x = letters[1:J],
                               cov7x = factor(1:J),
                               cov8x = c(TRUE, FALSE, TRUE)),
               repl_cov = NULL)
  expect_error(check_newdata(newdata, fit_no_repl_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov1x, cov2x, cov3x, cov4x \\(newdata\\); cov1, cov2, cov3, cov4 \\(fitted data\\)
  site_cov: cov5x, cov6x, cov7x, cov8x \\(newdata\\); cov5, cov6, cov7, cov8 \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov2 = rnorm(I),
                               cov1 = letters[1:I],
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov6 = rnorm(J),
                               cov5 = letters[1:J],
                               cov7 = factor(1:J),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = NULL)
  expect_error(check_newdata(newdata, fit_no_repl_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov2, cov1, cov3, cov4 \\(newdata\\); cov1, cov2, cov3, cov4 \\(fitted data\\)
  site_cov: cov6, cov5, cov7, cov8 \\(newdata\\); cov5, cov6, cov7, cov8 \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1x = rnorm(I),
                               cov2x = letters[1:I],
                               cov3x = factor(1:I),
                               cov4x = c(TRUE, FALSE)),
               site_cov = list(cov5x = rnorm(J),
                               cov6x = letters[1:J],
                               cov7x = factor(1:J),
                               cov8x = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K)))
  expect_error(check_newdata(newdata, fit_no_repl_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov1x, cov2x, cov3x, cov4x \\(newdata\\); cov1, cov2, cov3, cov4 \\(fitted data\\)
  site_cov: cov5x, cov6x, cov7x, cov8x \\(newdata\\); cov5, cov6, cov7, cov8 \\(fitted data\\)
  repl_cov: cov9 \\(newdata\\); \\(None\\) \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov2 = rnorm(I),
                               cov1 = letters[1:I],
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov6 = rnorm(J),
                               cov5 = letters[1:J],
                               cov7 = factor(1:J),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K)))
  expect_error(check_newdata(newdata, fit_no_repl_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov2, cov1, cov3, cov4 \\(newdata\\); cov1, cov2, cov3, cov4 \\(fitted data\\)
  site_cov: cov6, cov5, cov7, cov8 \\(newdata\\); cov5, cov6, cov7, cov8 \\(fitted data\\)
  repl_cov: cov9 \\(newdata\\); \\(None\\) \\(fitted data\\)")

  # Check for covariate classes
  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = letters[1:I],
                               cov2 = rnorm(I),
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J],
                               cov7 = c(TRUE, FALSE, TRUE),
                               cov8 = factor(1:J)),
               repl_cov = NULL)
  expect_error(check_newdata(newdata, fit_no_repl_cov),
               "The covariate classes in 'newdata' must match those in the fitted data.
  cov1: character \\(newdata\\), numeric \\(fitted data\\)
  cov2: numeric \\(newdata\\), character \\(fitted data\\)
  cov7: logical \\(newdata\\), factor \\(fitted data\\)
  cov8: factor \\(newdata\\), logical \\(fitted data\\)")

  # Check for new levels in discrete covariates
  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = rnorm(I),
                               cov2 = letters[1:I + 1],
                               cov3 = factor(1:I + 1),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J + 1],
                               cov7 = factor(1:J + 1),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = NULL)
  expect_error(check_newdata(newdata, fit_no_repl_cov),
               "The levels of discrete covariates in 'newdata' must match those in the fitted data.
  cov2: b, c \\(newdata\\); a, b \\(fitted data\\)
  cov3: 2, 3 \\(newdata\\); 1, 2 \\(fitted data\\)
  cov6: b, c, d \\(newdata\\); a, b, c \\(fitted data\\)
  cov7: 2, 3, 4 \\(newdata\\); 1, 2, 3 \\(fitted data\\)")

  # Check for orders of factor levels
  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = rnorm(I),
                               cov2 = letters[1:I],
                               cov3 = factor(1:I, levels = I:1),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J],
                               cov7 = factor(1:J, levels = J:1),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = NULL)
  expect_error(check_newdata(newdata, fit_no_repl_cov),
               "The levels of discrete covariates in 'newdata' must match those in the fitted data.
  cov3: 1, 2 \\(newdata\\); 1, 2 \\(fitted data\\)
  cov7: 1, 2, 3 \\(newdata\\); 1, 2, 3 \\(fitted data\\)")

  ## Extra covariates in newdata
  data_no_cov <- occumbData(y = y,
                            spec_cov = NULL,
                            site_cov = NULL,
                            repl_cov = NULL)

  fit_no_cov <-
    occumb(data = data_no_cov,
           n.chains = 1, n.burnin = 10, n.thin = 1, n.iter = 20,
           verbose = FALSE)

  expect_no_error(check_newdata(data_no_cov, fit_no_cov))

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = rnorm(I),
                               cov2 = letters[1:I],
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = NULL,
               repl_cov = NULL)
  expect_error(check_newdata(newdata, fit_no_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov1, cov2, cov3, cov4 \\(newdata\\); \\(None\\) \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = NULL,
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J],
                               cov7 = factor(1:J),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = NULL)
  expect_error(check_newdata(newdata, fit_no_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  site_cov: cov5, cov6, cov7, cov8 \\(newdata\\); \\(None\\) \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = NULL,
               site_cov = NULL,
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))

  expect_error(check_newdata(newdata, fit_no_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  repl_cov: cov9, cov10, cov11 \\(newdata\\); \\(None\\) \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = rnorm(I),
                               cov2 = letters[1:I],
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J],
                               cov7 = factor(1:J),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = NULL)
  expect_error(check_newdata(newdata, fit_no_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov1, cov2, cov3, cov4 \\(newdata\\); \\(None\\) \\(fitted data\\)
  site_cov: cov5, cov6, cov7, cov8 \\(newdata\\); \\(None\\) \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = rnorm(I),
                               cov2 = letters[1:I],
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = NULL,
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit_no_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov1, cov2, cov3, cov4 \\(newdata\\); \\(None\\) \\(fitted data\\)
  repl_cov: cov9, cov10, cov11 \\(newdata\\); \\(None\\) \\(fitted data\\)")

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = rnorm(I),
                               cov2 = letters[1:I],
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J],
                               cov7 = factor(1:J),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))
  expect_error(check_newdata(newdata, fit_no_cov),
               "The names of the covariates in 'newdata' and their order must match those in the fitted data.
  spec_cov: cov1, cov2, cov3, cov4 \\(newdata\\); \\(None\\) \\(fitted data\\)
  site_cov: cov5, cov6, cov7, cov8 \\(newdata\\); \\(None\\) \\(fitted data\\)
  repl_cov: cov9, cov10, cov11 \\(newdata\\); \\(None\\) \\(fitted data\\)")
})

### Tests for newdata labels ---------------------------------------------------
test_that("Species labels are copied correctly", {
  fit <- occumb(data = data,
                n.chains = 1, n.burnin = 10, n.thin = 1, n.iter = 20,
                verbose = FALSE)

  newdata <-
    occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)),
               spec_cov = list(cov1 = rnorm(I),
                               cov2 = letters[1:I],
                               cov3 = factor(1:I),
                               cov4 = c(TRUE, FALSE)),
               site_cov = list(cov5 = rnorm(J),
                               cov6 = letters[1:J],
                               cov7 = factor(1:J),
                               cov8 = c(TRUE, FALSE, TRUE)),
               repl_cov = list(cov9 = matrix(rnorm(J * K), J, K),
                               cov10 = matrix(letters[1:(J * K)], J, K),
                               cov11 = matrix(rep(c(TRUE, FALSE), J * K / 2), J, K)))

  expect_identical(attributes(suppressWarnings(predict(fit, newdata)))$label$Species,
                   dimnames(fit@data@y)[[1]])
})

### Tests for the function output ----------------------------------------------
test_that("Prediction and addition of attributes for phi works correctly", {
  N <- 10

  ## No shared effects, type = "i"
  fit <- occumb(data = data,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = N,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_phi),
                             formula(fit@occumb_args$formula_phi_shared),
                             "phi")
  post_effect <- get_post_samples(fit, "alpha")

  pred_link <- matrix(nrow = N, ncol = I)
  for (i in seq_len(I)) {
    pred_link[, i] <- post_effect[, i, ] * list_cov$cov
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "phi", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "phi", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- exp(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "phi", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "phi", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(exp(pred_link), 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "phi", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "phi", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- exp(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  ## No shared effects, type = "ij"
  fit <- occumb(data = data, formula_phi = ~ cov5,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_phi),
                             formula(fit@occumb_args$formula_phi_shared),
                             "phi")
  post_effect <- get_post_samples(fit, "alpha")

  pred_link <- array(dim = c(N, I, J))
  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      pred_link[, i, j] <- post_effect[, i, ] %*% list_cov$cov[j, ]
    }
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "phi", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "phi", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- exp(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "phi", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "phi", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- array(apply(exp(pred_link), 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "phi", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I, J))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "phi", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I, J))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- exp(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  ## No shared effects, type = "ijk"
  fit <- occumb(data = data, formula_phi = ~ cov9,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_phi),
                             formula(fit@occumb_args$formula_phi_shared),
                             "phi")
  post_effect <- get_post_samples(fit, "alpha")

  pred_link <- array(dim = c(N, I, J, K))
  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      for (k in seq_len(K)) {
        pred_link[, i, j, k] <- post_effect[, i, ] %*% list_cov$cov[j, k, ]
      }
    }
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "phi", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "phi", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- exp(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "phi", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "phi", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- array(apply(exp(pred_link), 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "phi", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Samples", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "phi", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Samples", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- exp(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  ## With shared effects, type = "i"
  fit <- occumb(data = data, formula_phi_shared = ~ cov1,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_phi),
                             formula(fit@occumb_args$formula_phi_shared),
                             "phi")
  post_effect <- get_post_samples(fit, "alpha")
  post_effect_shared <- get_post_samples(fit, "alpha_shared")

  pred_link <- matrix(nrow = N, ncol = I)
  for (i in seq_len(I)) {
    pred_link[, i] <-
      post_effect[, i, ] * list_cov$cov +
      post_effect_shared * list_cov$cov_shared[i, ]
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "phi", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "phi", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- exp(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "phi", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "phi", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(exp(pred_link), 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "phi", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "phi", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- exp(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  fit <- occumb(data = data, formula_phi_shared = ~ cov1 + cov2,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_phi),
                             formula(fit@occumb_args$formula_phi_shared),
                             "phi")
  post_effect <- get_post_samples(fit, "alpha")
  post_effect_shared <- get_post_samples(fit, "alpha_shared")

  pred_link <- matrix(nrow = N, ncol = I)
  for (i in seq_len(I)) {
    pred_link[, i] <-
      post_effect[, i, ] * list_cov$cov +
      post_effect_shared %*% list_cov$cov_shared[i, ]
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "phi", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "phi", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- exp(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "phi", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "phi", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(exp(pred_link), 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "phi", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "phi", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- exp(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  ## With shared effects, type = "ij"
  fit <- occumb(data = data, formula_phi = ~ cov5, formula_phi_shared = ~ cov1,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_phi),
                             formula(fit@occumb_args$formula_phi_shared),
                             "phi")
  post_effect <- get_post_samples(fit, "alpha")
  post_effect_shared <- get_post_samples(fit, "alpha_shared")

  pred_link <- array(dim = c(N, I, J))
  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      pred_link[, i, j] <-
        post_effect[, i, ] %*% list_cov$cov[j, ] +
        post_effect_shared * list_cov$cov_shared[i, j, ]
    }
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "phi", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "phi", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- exp(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "phi", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "phi", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- array(apply(exp(pred_link), 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "phi", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I, J))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "phi", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I, J))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- exp(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  ## With shared effects, type = "ijk"
  fit <- occumb(data = data, formula_phi = ~ cov9, formula_phi_shared = ~ cov1,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_phi),
                             formula(fit@occumb_args$formula_phi_shared),
                             "phi")
  post_effect <- get_post_samples(fit, "alpha")
  post_effect_shared <- get_post_samples(fit, "alpha_shared")

  pred_link <- array(dim = c(N, I, J, K))
  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      for (k in seq_len(K)) {
        pred_link[, i, j, k] <-
          post_effect[, i, ] %*% list_cov$cov[j, k, ] +
          post_effect_shared * list_cov$cov_shared[i, j, k, ]
      }
    }
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "phi", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "phi", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- exp(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "phi", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "phi", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- array(apply(exp(pred_link), 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "phi", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Samples", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "phi", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Samples", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- exp(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)
})

test_that("Prediction and addition of attributes for theta works correctly", {
  N <- 10

  ## No shared effects, type = "i"
  fit <- occumb(data = data,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = N,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_theta),
                             formula(fit@occumb_args$formula_theta_shared),
                             "theta")
  post_effect <- get_post_samples(fit, "beta")

  pred_link <- matrix(nrow = N, ncol = I)
  for (i in seq_len(I)) {
    pred_link[, i] <- post_effect[, i, ] * list_cov$cov
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "theta", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "theta", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- plogis(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "theta", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "theta", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(plogis(pred_link), 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "theta", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "theta", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- plogis(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  ## No shared effects, type = "ij"
  fit <- occumb(data = data, formula_theta = ~ cov5,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_theta),
                             formula(fit@occumb_args$formula_theta_shared),
                             "theta")
  post_effect <- get_post_samples(fit, "beta")

  pred_link <- array(dim = c(N, I, J))
  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      pred_link[, i, j] <- post_effect[, i, ] %*% list_cov$cov[j, ]
    }
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "theta", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "theta", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- plogis(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "theta", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "theta", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- array(apply(plogis(pred_link), 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "theta", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I, J))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "theta", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I, J))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- plogis(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  ## No shared effects, type = "ijk"
  fit <- occumb(data = data, formula_theta = ~ cov9,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_theta),
                             formula(fit@occumb_args$formula_theta_shared),
                             "theta")
  post_effect <- get_post_samples(fit, "beta")

  pred_link <- array(dim = c(N, I, J, K))
  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      for (k in seq_len(K)) {
        pred_link[, i, j, k] <- post_effect[, i, ] %*% list_cov$cov[j, k, ]
      }
    }
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "theta", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "theta", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- plogis(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "theta", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "theta", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- array(apply(plogis(pred_link), 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "theta", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Samples", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "theta", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Samples", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- plogis(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  ## With shared effects, type = "i"
  fit <- occumb(data = data, formula_theta_shared = ~ cov1,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_theta),
                             formula(fit@occumb_args$formula_theta_shared),
                             "theta")
  post_effect <- get_post_samples(fit, "beta")
  post_effect_shared <- get_post_samples(fit, "beta_shared")

  pred_link <- matrix(nrow = N, ncol = I)
  for (i in seq_len(I)) {
    pred_link[, i] <-
      post_effect[, i, ] * list_cov$cov +
      post_effect_shared * list_cov$cov_shared[i, ]
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "theta", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "theta", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- plogis(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "theta", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "theta", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(plogis(pred_link), 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "theta", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "theta", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- plogis(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  fit <- occumb(data = data, formula_theta_shared = ~ cov1 + cov2,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_theta),
                             formula(fit@occumb_args$formula_theta_shared),
                             "theta")
  post_effect <- get_post_samples(fit, "beta")
  post_effect_shared <- get_post_samples(fit, "beta_shared")

  pred_link <- matrix(nrow = N, ncol = I)
  for (i in seq_len(I)) {
    pred_link[, i] <-
      post_effect[, i, ] * list_cov$cov +
      post_effect_shared %*% list_cov$cov_shared[i, ]
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "theta", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "theta", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- plogis(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "theta", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "theta", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(plogis(pred_link), 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "theta", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "theta", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- plogis(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  ## With shared effects, type = "ij"
  fit <- occumb(data = data, formula_theta = ~ cov5, formula_theta_shared = ~ cov1,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_theta),
                             formula(fit@occumb_args$formula_theta_shared),
                             "theta")
  post_effect <- get_post_samples(fit, "beta")
  post_effect_shared <- get_post_samples(fit, "beta_shared")

  pred_link <- array(dim = c(N, I, J))
  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      pred_link[, i, j] <-
        post_effect[, i, ] %*% list_cov$cov[j, ] +
        post_effect_shared * list_cov$cov_shared[i, j, ]
    }
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "theta", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "theta", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- plogis(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "theta", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "theta", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- array(apply(plogis(pred_link), 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "theta", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I, J))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "theta", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I, J))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- plogis(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  ## With shared effects, type = "ijk"
  fit <- occumb(data = data, formula_theta = ~ cov9, formula_theta_shared = ~ cov1,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_theta),
                             formula(fit@occumb_args$formula_theta_shared),
                             "theta")
  post_effect <- get_post_samples(fit, "beta")
  post_effect_shared <- get_post_samples(fit, "beta_shared")

  pred_link <- array(dim = c(N, I, J, K))
  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      for (k in seq_len(K)) {
        pred_link[, i, j, k] <-
          post_effect[, i, ] %*% list_cov$cov[j, k, ] +
          post_effect_shared * list_cov$cov_shared[i, j, k, ]
      }
    }
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "theta", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "theta", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- plogis(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "theta", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "theta", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Statistics", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- array(apply(plogis(pred_link), 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "theta", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Samples", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "theta", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I, J, K))
  expect_equal(attributes(pred)$dimension,
               c("Samples", "Species", "Sites", "Replicates"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])
  expect_equal(attributes(pred)$label$Replicates, dimnames(fit@data@y)[[3]])

  ans  <- plogis(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)
})

test_that("Prediction and addition of attributes for psi works correctly", {
  N <- 10

  ## No shared effects, type = "i"
  fit <- occumb(data = data,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = N,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_psi),
                             formula(fit@occumb_args$formula_psi_shared),
                             "psi")
  post_effect <- get_post_samples(fit, "gamma")

  pred_link <- matrix(nrow = N, ncol = I)
  for (i in seq_len(I)) {
    pred_link[, i] <- post_effect[, i, ] * list_cov$cov
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "psi", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "psi", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- plogis(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "psi", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "psi", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(plogis(pred_link), 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "psi", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "psi", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- plogis(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  ## No shared effects, type = "ij"
  fit <- occumb(data = data, formula_psi = ~ cov5,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_psi),
                             formula(fit@occumb_args$formula_psi_shared),
                             "psi")
  post_effect <- get_post_samples(fit, "gamma")

  pred_link <- array(dim = c(N, I, J))
  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      pred_link[, i, j] <- post_effect[, i, ] %*% list_cov$cov[j, ]
    }
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "psi", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "psi", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- plogis(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "psi", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "psi", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- array(apply(plogis(pred_link), 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "psi", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I, J))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "psi", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I, J))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- plogis(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  ## With shared effects, type = "i"
  fit <- occumb(data = data, formula_psi_shared = ~ cov1,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_psi),
                             formula(fit@occumb_args$formula_psi_shared),
                             "psi")
  post_effect <- get_post_samples(fit, "gamma")
  post_effect_shared <- get_post_samples(fit, "gamma_shared")

  pred_link <- matrix(nrow = N, ncol = I)
  for (i in seq_len(I)) {
    pred_link[, i] <-
      post_effect[, i, ] * list_cov$cov +
      post_effect_shared * list_cov$cov_shared[i, ]
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "psi", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "psi", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- plogis(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "psi", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "psi", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(plogis(pred_link), 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "psi", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "psi", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- plogis(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  fit <- occumb(data = data, formula_psi_shared = ~ cov1 + cov2,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_psi),
                             formula(fit@occumb_args$formula_psi_shared),
                             "psi")
  post_effect <- get_post_samples(fit, "gamma")
  post_effect_shared <- get_post_samples(fit, "gamma_shared")

  pred_link <- matrix(nrow = N, ncol = I)
  for (i in seq_len(I)) {
    pred_link[, i] <-
      post_effect[, i, ] * list_cov$cov +
      post_effect_shared %*% list_cov$cov_shared[i, ]
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "psi", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "psi", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- plogis(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "psi", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "psi", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- array(apply(plogis(pred_link), 2:length(dim(pred_link)), mean), c(1, ncol(pred_link)))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "psi", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "psi", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])

  ans  <- plogis(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)


  ## With shared effects, type = "ij"
  fit <- occumb(data = data, formula_psi = ~ cov5, formula_psi_shared = ~ cov1,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = 10,
                verbose = FALSE)

  list_cov <- set_covariates(data,
                             formula(fit@occumb_args$formula_psi),
                             formula(fit@occumb_args$formula_psi_shared),
                             "psi")
  post_effect <- get_post_samples(fit, "gamma")
  post_effect_shared <- get_post_samples(fit, "gamma_shared")

  pred_link <- array(dim = c(N, I, J))
  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      pred_link[, i, j] <-
        post_effect[, i, ] %*% list_cov$cov[j, ] +
        post_effect_shared * list_cov$cov_shared[i, j, ]
    }
  }

  # 1. scale = link, type = quantiles
  pred <- predict(fit, parameter = "psi", scale = "link", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- apply(pred_link, 2:length(dim(pred_link)),
                quantile, probs = c(0.5, 0.025, 0.975))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 2. scale = response, type = quantiles
  pred <- predict(fit, parameter = "psi", scale = "response", type = "quantiles")

  expect_equal(dim(pred), c(3, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, c("50%", "2.5%", "97.5%"))
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- plogis(ans)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 3. scale = link, type = mean
  pred <- predict(fit, parameter = "psi", scale = "link", type = "mean")

  expect_equal(dim(pred), c(1, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- array(apply(pred_link, 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 4. scale = response, type = mean
  pred <- predict(fit, parameter = "psi", scale = "response", type = "mean")

  expect_equal(dim(pred), c(1, I, J))
  expect_equal(attributes(pred)$dimension, c("Statistics", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Statistics, "mean")
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- array(apply(plogis(pred_link), 2:length(dim(pred_link)), mean),
                c(1, dim(pred_link)[-1]))
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 5. scale = link, type = samples
  pred <- predict(fit, parameter = "psi", scale = "link", type = "samples")

  expect_equal(dim(pred), c(N, I, J))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- pred_link
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)

  # 6. scale = response, type = samples
  pred <- predict(fit, parameter = "psi", scale = "response", type = "samples")

  expect_equal(dim(pred), c(N, I, J))
  expect_equal(attributes(pred)$dimension, c("Samples", "Species", "Sites"))
  expect_equal(attributes(pred)$label$Samples, NULL)
  expect_equal(attributes(pred)$label$Species, dimnames(fit@data@y)[[1]])
  expect_equal(attributes(pred)$label$Sites, dimnames(fit@data@y)[[2]])

  ans  <- plogis(pred_link)
  attr(pred, "parameter") <- attr(pred, "scale") <-
    attr(pred, "dimension") <- attr(pred, "label") <- NULL

  expect_equal(pred, ans)
})
