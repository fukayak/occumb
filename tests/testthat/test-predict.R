### Test data ------------------------------------------------------------------
I <- 2
J <- 3
K <- 4
N <- 5
y <- y_unnamed <- array(sample.int(I * J * K), dim = c(I, J, K))
dimnames(y)[[1]] <- sprintf("species %s", 1:I)
dimnames(y)[[2]] <- sprintf("site %s", 1:J)
dimnames(y)[[3]] <- sprintf("replicate %s", 1:K)
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
data_unnamed <-
  occumbData(y = y_unnamed,
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
test_that("Prediction and addition of attributes works correctly", {
  test_array_output <- function(data,
                                formula = ~ 1,
                                formula_shared = ~ 1,
                                test_psi = TRUE,
                                n.chains = 1,
                                n.burnin = 0,
                                n.thin = 1,
                                n.iter) {

    if (test_psi) {
      fit <- occumb(data = data,
                    formula_phi = formula,
                    formula_theta = formula,
                    formula_psi = formula,
                    formula_phi_shared = formula_shared,
                    formula_theta_shared = formula_shared,
                    formula_psi_shared = formula_shared,
                    n.chains = n.chains,
                    n.burnin = n.burnin,
                    n.thin = n.thin,
                    n.iter = n.iter,
                    verbose = FALSE)
    } else {
      fit <- occumb(data = data,
                    formula_phi = formula,
                    formula_theta = formula,
                    formula_psi = ~ 1,
                    formula_phi_shared = formula_shared,
                    formula_theta_shared = formula_shared,
                    formula_psi_shared = ~ 1,
                    n.chains = n.chains,
                    n.burnin = n.burnin,
                    n.thin = n.thin,
                    n.iter = n.iter,
                    verbose = FALSE)
    }

    .test_array_output(fit, "phi", "link",     "quantiles")
    .test_array_output(fit, "phi", "link",     "mean")
    .test_array_output(fit, "phi", "link",     "samples")
    .test_array_output(fit, "phi", "response", "quantiles")
    .test_array_output(fit, "phi", "response", "mean")
    .test_array_output(fit, "phi", "response", "samples")

    .test_array_output(fit, "theta", "link",     "quantiles")
    .test_array_output(fit, "theta", "link",     "mean")
    .test_array_output(fit, "theta", "link",     "samples")
    .test_array_output(fit, "theta", "response", "quantiles")
    .test_array_output(fit, "theta", "response", "mean")
    .test_array_output(fit, "theta", "response", "samples")

    if (test_psi) {
      .test_array_output(fit, "psi", "link",     "quantiles")
      .test_array_output(fit, "psi", "link",     "mean")
      .test_array_output(fit, "psi", "link",     "samples")
      .test_array_output(fit, "psi", "response", "quantiles")
      .test_array_output(fit, "psi", "response", "mean")
      .test_array_output(fit, "psi", "response", "samples")
    }
  }

  .test_array_output <- function(fit, parameter, scale, type) {

    test_object_structure <- function(fit, parameter, scale, type) {

      pred <- predict(fit, parameter = parameter, scale = scale, type = type)
      parameter_dimension <- set_parameter_dimension(fit, parameter)
      dim_ans <- set_dim_ans(type, parameter_dimension)

      expect_equal(dim(pred), dim_ans)
    }

    test_object_attributes <- function(fit, parameter, scale, type) {

      pred <- predict(fit, parameter = parameter, scale = scale, type = type)
      parameter_dimension <- set_parameter_dimension(fit, parameter)
      label_dimension_ans <-
        set_label_dimension_ans(type, parameter_dimension)
      label_statistics_ans <- set_label_statistics_ans(type)

      expect_equal(attributes(pred)$dimension, label_dimension_ans)

      if (type == "samples") {
        expect_equal(attributes(pred)$label$Samples, NULL)
      } else {
        expect_equal(attributes(pred)$label$Statistics, label_statistics_ans)
      }

      expect_equal(attributes(pred)$label$Species,
                   dimnames(fit@data@y)[[1]])
      if (parameter_dimension != "i") {
        expect_equal(attributes(pred)$label$Sites,
                     dimnames(fit@data@y)[[2]])
      }
      if (parameter_dimension == "ijk") {
        expect_equal(attributes(pred)$label$Replicates,
                     dimnames(fit@data@y)[[3]])
      }
    }

    test_object_contents <- function(fit, parameter, scale, type) {

      pred <- predict(fit, parameter = parameter, scale = scale, type = type)
      pred_without_attributes <- pred |>
        remove_attributes("parameter") |>
        remove_attributes("scale") |>
        remove_attributes("dimension") |>
        remove_attributes("label")
      pred_ans <- set_pred_ans(fit, parameter, scale, type)

      expect_equal(pred_without_attributes, pred_ans)
    }

    set_parameter_dimension <- function(fit, parameter) {
      eval(parse(text = paste0("get_modargs(fit)$", parameter)))
    }

    set_dim_ans <- function(type, parameter_dimension) {
      if (type == "quantiles") {
        if (parameter_dimension == "i") {
          return(c(3, I))
        }
        if (parameter_dimension == "ij") {
          return(c(3, I, J))
        }
        if (parameter_dimension == "ijk") {
          return(c(3, I, J, K))
        }
      }

      if (type == "mean") {
        if (parameter_dimension == "i") {
          return(c(1, I))
        }
        if (parameter_dimension == "ij") {
          return(c(1, I, J))
        }
        if (parameter_dimension == "ijk") {
          return(c(1, I, J, K))
        }
      }

      if (type == "samples") {
        if (parameter_dimension == "i") {
          return(c(N, I))
        }
        if (parameter_dimension == "ij") {
          return(c(N, I, J))
        }
        if (parameter_dimension == "ijk") {
          return(c(N, I, J, K))
        }
      }
    }

    set_label_dimension_ans <- function(type, parameter_dimension) {
      if (type == "samples") {
        if (parameter_dimension == "i") {
          return(c("Samples", "Species"))
        }
        if (parameter_dimension == "ij") {
          return(c("Samples", "Species", "Sites"))
        }
        if (parameter_dimension == "ijk") {
          return(c("Samples", "Species", "Sites", "Replicates"))
        }
      } else {
        if (parameter_dimension == "i") {
          return(c("Statistics", "Species"))
        }
        if (parameter_dimension == "ij") {
          return(c("Statistics", "Species", "Sites"))
        }
        if (parameter_dimension == "ijk") {
          return(c("Statistics", "Species", "Sites", "Replicates"))
        }
      }
    }

    set_label_statistics_ans <- function(type) {
      if (type == "quantiles") {
        return(c("50%", "2.5%", "97.5%"))
      }
      if (type == "mean") {
        return("mean")
      }
    }

    remove_attributes <- function(x, attribute) {
      attr(x, attribute) <- NULL
      return(x)
    }

    set_pred_ans <- function(fit, parameter, scale, type) {
      pred_ans <- get_pred_link(fit, parameter, scale, type) |>
        get_pred_ans(fit, parameter, scale, type)

      return(pred_ans)
    }

    get_pred_link <- function(fit, parameter, scale, type) {

      parameter_dimension <- set_parameter_dimension(fit, parameter)
      has_shared_effect   <- set_has_shared_effect(fit, parameter)
      coefficient         <- set_coefficient(parameter)
      coefficient_shared  <- set_coefficient_shared(parameter)
      list_cov            <- get_covariates(fit, parameter)
      post_effect         <- get_post_samples(fit, coefficient)
      if (has_shared_effect) {
        post_effect_shared <- get_post_samples(fit, coefficient_shared)
      }

      if (parameter_dimension == "i") {
        pred_link <- matrix(nrow = N, ncol = I)

        if (!has_shared_effect) {
          for (i in seq_len(I)) {
            pred_link[, i] <- post_effect[, i, ] * list_cov$cov
          }
        } else {
          if (length(post_effect_shared) > 1) {
            for (i in seq_len(I)) {
              pred_link[, i] <- post_effect[, i, ] * list_cov$cov +
                post_effect_shared %*% list_cov$cov_shared[i, ]
            }
          } else {
            for (i in seq_len(I)) {
              pred_link[, i] <- post_effect[, i, ] * list_cov$cov +
                post_effect_shared * list_cov$cov_shared[i, ]
            }
          }
        }

        return(pred_link)
      }

      if (parameter_dimension == "ij") {
        pred_link <- array(dim = c(N, I, J))

        if (!has_shared_effect) {
          for (i in seq_len(I)) {
            for (j in seq_len(J)) {
              pred_link[, i, j] <- post_effect[, i, ] %*% list_cov$cov[j, ]
            }
          }
        } else {
          if (length(post_effect_shared) > 1) {
            for (i in seq_len(I)) {
              for (j in seq_len(J)) {
                pred_link[, i, j] <- post_effect[, i, ] %*% list_cov$cov[j, ] +
                  post_effect_shared %*% list_cov$cov_shared[i, j, ]
              }
            }
          } else {
            for (i in seq_len(I)) {
              for (j in seq_len(J)) {
                pred_link[, i, j] <- post_effect[, i, ] %*% list_cov$cov[j, ] +
                  post_effect_shared * list_cov$cov_shared[i, j, ]
              }
            }
          }
        }

        return(pred_link)
      }

      if (parameter_dimension == "ijk") {
        pred_link <- array(dim = c(N, I, J, K))

        if (!has_shared_effect) {
          for (i in seq_len(I)) {
            for (j in seq_len(J)) {
              for (k in seq_len(K)) {
                pred_link[, i, j, k] <- post_effect[, i, ] %*% list_cov$cov[j, k, ]
              }
            }
          }
        } else {
          if (length(post_effect_shared) > 1) {
            for (i in seq_len(I)) {
              for (j in seq_len(J)) {
                for (k in seq_len(K)) {
                  pred_link[, i, j, k] <-
                    post_effect[, i, ] %*% list_cov$cov[j, k, ] +
                    post_effect_shared %*% list_cov$cov_shared[i, j, k, ]
                }
              }
            }
          } else {
            for (i in seq_len(I)) {
              for (j in seq_len(J)) {
                for (k in seq_len(K)) {
                  pred_link[, i, j, k] <-
                    post_effect[, i, ] %*% list_cov$cov[j, k, ] +
                    post_effect_shared * list_cov$cov_shared[i, j, k, ]
                }
              }
            }
          }
        }

        return(pred_link)
      }
    }

    get_pred_ans <- function(pred_link, fit, parameter, scale, type) {

      parameter_dimension <- set_parameter_dimension(fit, parameter)
      inv_link <- set_inv_link(parameter)

      if (type == "quantiles") {
        pred_ans <- apply(pred_link, 2:length(dim(pred_link)),
                          quantile, probs = c(0.5, 0.025, 0.975))
        if (scale == "response") {
          pred_ans <- inv_link(pred_ans)
        }

        return(pred_ans)
      }

      if (type == "mean") {
        if (scale == "link") {
          if (parameter_dimension == "i") {
            pred_ans <- array(apply(pred_link, 2:length(dim(pred_link)), mean),
                              c(1, ncol(pred_link)))
          } else {
            pred_ans <- array(apply(pred_link, 2:length(dim(pred_link)), mean),
                              c(1, dim(pred_link)[-1]))
          }
        }
        if (scale == "response") {
          if (parameter_dimension == "i") {
            pred_ans <-
              array(apply(inv_link(pred_link), 2:length(dim(pred_link)), mean),
                    c(1, ncol(pred_link)))
          } else {
            pred_ans <-
              array(apply(inv_link(pred_link), 2:length(dim(pred_link)), mean),
                    c(1, dim(pred_link)[-1]))
          }
        }

        return(pred_ans)
      }

      if (type == "samples") {
        if (scale == "link") {
          pred_ans <- pred_link
        }
        if (scale == "response") {
          pred_ans <- inv_link(pred_link)
        }

        return(pred_ans)
      }
    }

    set_has_shared_effect <- function(fit, parameter) {
      eval(parse(text = paste0("get_modargs(fit)$", parameter, "_shared")))
    }

    set_coefficient <- function(parameter) {
      switch(parameter,
             phi = "alpha",
             theta = "beta",
             psi = "gamma")
    }

    set_coefficient_shared <- function(parameter) {
      switch(parameter,
             phi = "alpha_shared",
             theta = "beta_shared",
             psi = "gamma_shared")
    }

    set_inv_link <- function(parameter) {
      switch(parameter,
             phi   = exp,
             theta = stats::plogis,
             psi   = stats::plogis)
    }

    test_object_structure(fit, parameter, scale, type)
    test_object_attributes(fit, parameter, scale, type)
    test_object_contents(fit, parameter, scale, type)
  }

  ## No shared effects, type = "i"
  test_array_output(data, n.iter = N)

  ## No shared effects, type = "ij"
  test_array_output(data, formula = ~ cov5, n.iter = N)

  ## No shared effects, type = "ijk"
  test_array_output(data, formula = ~ cov9, test_psi = FALSE, n.iter = N)

  ## With shared effects, type = "i"
  test_array_output(data, formula_shared = ~ cov1, n.iter = N)

  # With discrete covariates
  test_array_output(data, formula_shared = ~ cov1 + cov2, n.iter = N)

  ## With shared effects, type = "ij"
  test_array_output(data, formula = ~ cov5, formula_shared = ~ cov1, n.iter = N)

  ## With shared effects, type = "ijk"
  test_array_output(data, formula = ~ cov9, formula_shared = ~ cov1, test_psi = FALSE,
      n.iter = N)
})


test_that("Option for dataframe output works", {

  test_dataframe_output <- function(fit, parameter, scale, type) {

    pred <- predict(fit, parameter = parameter, scale = scale, type = type,
                    output_dataframe = TRUE)
    pred_array <- predict(fit, parameter = parameter, scale = scale, type = type,
                          output_dataframe = FALSE)
    parameter_dimension <- 
      eval(parse(text = paste0("get_modargs(fit)$", parameter)))

    multiplier_type <- switch(type,
                              quantiles = 3,
                              mean = 1,
                              samples = N)
    multiplier_parameter_dimension <- switch(parameter_dimension,
                                             i = I,
                                             ij = I * J,
                                             ijk = I * J * K)
    colnames_dimension <- switch(parameter_dimension,
                                 i = c("Species"),
                                 ij = c("Species", "Sites"),
                                 ijk = c("Species", "Sites", "Replicates"))

    nrows_ans <- multiplier_type * multiplier_parameter_dimension
    ncols_ans <- switch(parameter_dimension,
                        i = 5,
                        ij = 6,
                        ijk = 7)
    colnames_ans <- if (type == "samples") {
      c(c("Parameter", "Scale", "Samples"), colnames_dimension, "Value")
    } else {
      c(c("Parameter", "Scale", "Statistics"), colnames_dimension, "Value")
    }

    Parameter_ans <- factor(rep(parameter, nrows_ans))
    Scale_ans <- factor(rep(attributes(pred_array)$scale, nrows_ans))
    Statistics_Samples_ans <- if (type == "samples") {
      factor(rep(seq_len(N), multiplier_parameter_dimension))
    } else {
      factor(rep(attributes(pred_array)$label$Statistics,
                 multiplier_parameter_dimension),
             levels = attributes(pred_array)$label$Statistics)
    }

    # Tests for object structure
    checkmate::expect_data_frame(pred, nrows = nrows_ans, ncols = ncols_ans)
    expect_identical(colnames(pred), colnames_ans)

    # Tests for object contents
    expect_identical(pred$Parameter, Parameter_ans)
    expect_identical(pred$Scale, Scale_ans)
    expect_identical(pred[[3]], Statistics_Samples_ans)
    expect_identical(pred$Value, c(pred_array))

    if (is.null(attributes(pred_array)$label$Species)) {
      label_spec <- seq_len(I)
    } else {
      label_spec <- attributes(pred_array)$label$Species
    }

    if (is.null(attributes(pred_array)$label$Sites)) {
      label_site <- seq_len(J)
    } else {
      label_site <- attributes(pred_array)$label$Sites
    }

    if (is.null(attributes(pred_array)$label$Replicates)) {
      label_repl <- seq_len(K)
    } else {
      label_repl <- attributes(pred_array)$label$Replicates
    }

    if (parameter_dimension == "i") {
      labels_ans <- expand.grid(seq_len(multiplier_type),
                                label_spec)

      Species_ans <- factor(labels_ans[, 2])
      if (!is.null(attributes(pred_array)$label$Species))
        levels(Species_ans) <- attributes(pred_array)$label$Species

      expect_identical(pred$Species, Species_ans)
    }
    if (parameter_dimension == "ij") {
      labels_ans <- expand.grid(seq_len(multiplier_type),
                                label_spec,
                                label_site)

      Species_ans <- factor(labels_ans[, 2])
      if (!is.null(attributes(pred_array)$label$Species))
        levels(Species_ans) <- attributes(pred_array)$label$Species

      expect_identical(pred$Species, Species_ans)

      Sites_ans <- factor(labels_ans[, 3])
      if (!is.null(attributes(pred_array)$label$Sites))
        levels(Sites_ans) <- attributes(pred_array)$label$Sites

      expect_identical(pred$Sites, Sites_ans)
    }
    if (parameter_dimension == "ijk") {
      labels_ans <- expand.grid(seq_len(multiplier_type),
                                label_spec,
                                label_site,
                                label_repl)

      Species_ans <- factor(labels_ans[, 2])
      if (!is.null(attributes(pred_array)$label$Species))
        levels(Species_ans) <- attributes(pred_array)$label$Species

      expect_identical(pred$Species, Species_ans)

      Sites_ans <- factor(labels_ans[, 3])
      if (!is.null(attributes(pred_array)$label$Sites))
        levels(Sites_ans) <- attributes(pred_array)$label$Sites

      expect_identical(pred$Sites, Sites_ans)

      Replicates_ans <- factor(labels_ans[, 4])
      if (!is.null(attributes(pred_array)$label$Replicates))
        levels(Replicates_ans) <- attributes(pred_array)$label$Replicates

      expect_identical(pred$Replicates, Replicates_ans)
    }
  }

  ## type = "i"
  fit <- occumb(data = data,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = N,
                verbose = FALSE)

  fit_unnamed <- occumb(data = data_unnamed,
                        n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = N,
                        verbose = FALSE)

  # Tests for phi
  test_dataframe_output(fit,         "phi", "link", "quantiles")
  test_dataframe_output(fit,         "phi", "link", "mean")
  test_dataframe_output(fit,         "phi", "link", "samples")
  test_dataframe_output(fit,         "phi", "response", "quantiles")
  test_dataframe_output(fit,         "phi", "response", "mean")
  test_dataframe_output(fit,         "phi", "response", "samples")
  test_dataframe_output(fit_unnamed, "phi", "link", "quantiles")
  test_dataframe_output(fit_unnamed, "phi", "link", "mean")
  test_dataframe_output(fit_unnamed, "phi", "link", "samples")
  test_dataframe_output(fit_unnamed, "phi", "response", "quantiles")
  test_dataframe_output(fit_unnamed, "phi", "response", "mean")
  test_dataframe_output(fit_unnamed, "phi", "response", "samples")

  # Tests for theta
  test_dataframe_output(fit,         "theta", "link", "quantiles")
  test_dataframe_output(fit,         "theta", "link", "mean")
  test_dataframe_output(fit,         "theta", "link", "samples")
  test_dataframe_output(fit,         "theta", "response", "quantiles")
  test_dataframe_output(fit,         "theta", "response", "mean")
  test_dataframe_output(fit,         "theta", "response", "samples")
  test_dataframe_output(fit_unnamed, "theta", "link", "quantiles")
  test_dataframe_output(fit_unnamed, "theta", "link", "mean")
  test_dataframe_output(fit_unnamed, "theta", "link", "samples")
  test_dataframe_output(fit_unnamed, "theta", "response", "quantiles")
  test_dataframe_output(fit_unnamed, "theta", "response", "mean")
  test_dataframe_output(fit_unnamed, "theta", "response", "samples")

  # Tests for psi
  test_dataframe_output(fit,         "psi", "link", "quantiles")
  test_dataframe_output(fit,         "psi", "link", "mean")
  test_dataframe_output(fit,         "psi", "link", "samples")
  test_dataframe_output(fit,         "psi", "response", "quantiles")
  test_dataframe_output(fit,         "psi", "response", "mean")
  test_dataframe_output(fit,         "psi", "response", "samples")
  test_dataframe_output(fit_unnamed, "psi", "link", "quantiles")
  test_dataframe_output(fit_unnamed, "psi", "link", "mean")
  test_dataframe_output(fit_unnamed, "psi", "link", "samples")
  test_dataframe_output(fit_unnamed, "psi", "response", "quantiles")
  test_dataframe_output(fit_unnamed, "psi", "response", "mean")
  test_dataframe_output(fit_unnamed, "psi", "response", "samples")

  ## type = "ij"
  fit <- occumb(data = data,
                formula_phi = ~ cov5,
                formula_theta = ~ cov5,
                formula_psi = ~ cov5,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = N,
                verbose = FALSE)

  fit_unnamed <- occumb(data = data_unnamed,
                        formula_phi = ~ cov5,
                        formula_theta = ~ cov5,
                        formula_psi = ~ cov5,
                        n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = N,
                        verbose = FALSE)

  # Tests for phi
  test_dataframe_output(fit,         "phi", "link", "quantiles")
  test_dataframe_output(fit,         "phi", "link", "mean")
  test_dataframe_output(fit,         "phi", "link", "samples")
  test_dataframe_output(fit,         "phi", "response", "quantiles")
  test_dataframe_output(fit,         "phi", "response", "mean")
  test_dataframe_output(fit,         "phi", "response", "samples")
  test_dataframe_output(fit_unnamed, "phi", "link", "quantiles")
  test_dataframe_output(fit_unnamed, "phi", "link", "mean")
  test_dataframe_output(fit_unnamed, "phi", "link", "samples")
  test_dataframe_output(fit_unnamed, "phi", "response", "quantiles")
  test_dataframe_output(fit_unnamed, "phi", "response", "mean")
  test_dataframe_output(fit_unnamed, "phi", "response", "samples")

  # Tests for theta
  test_dataframe_output(fit,         "theta", "link", "quantiles")
  test_dataframe_output(fit,         "theta", "link", "mean")
  test_dataframe_output(fit,         "theta", "link", "samples")
  test_dataframe_output(fit,         "theta", "response", "quantiles")
  test_dataframe_output(fit,         "theta", "response", "mean")
  test_dataframe_output(fit,         "theta", "response", "samples")
  test_dataframe_output(fit_unnamed, "theta", "link", "quantiles")
  test_dataframe_output(fit_unnamed, "theta", "link", "mean")
  test_dataframe_output(fit_unnamed, "theta", "link", "samples")
  test_dataframe_output(fit_unnamed, "theta", "response", "quantiles")
  test_dataframe_output(fit_unnamed, "theta", "response", "mean")
  test_dataframe_output(fit_unnamed, "theta", "response", "samples")

  # Tests for psi
  test_dataframe_output(fit,         "psi", "link", "quantiles")
  test_dataframe_output(fit,         "psi", "link", "mean")
  test_dataframe_output(fit,         "psi", "link", "samples")
  test_dataframe_output(fit,         "psi", "response", "quantiles")
  test_dataframe_output(fit,         "psi", "response", "mean")
  test_dataframe_output(fit,         "psi", "response", "samples")
  test_dataframe_output(fit_unnamed, "psi", "link", "quantiles")
  test_dataframe_output(fit_unnamed, "psi", "link", "mean")
  test_dataframe_output(fit_unnamed, "psi", "link", "samples")
  test_dataframe_output(fit_unnamed, "psi", "response", "quantiles")
  test_dataframe_output(fit_unnamed, "psi", "response", "mean")
  test_dataframe_output(fit_unnamed, "psi", "response", "samples")

  ## type = "ijk"
  fit <- occumb(data = data,
                formula_phi = ~ cov9,
                formula_theta = ~ cov9,
                n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = N,
                verbose = FALSE)

  fit_unnamed <- occumb(data = data_unnamed,
                        formula_phi = ~ cov9,
                        formula_theta = ~ cov9,
                        n.chains = 1, n.burnin = 0, n.thin = 1, n.iter = N,
                        verbose = FALSE)

  # Tests for phi
  test_dataframe_output(fit,         "phi", "link", "quantiles")
  test_dataframe_output(fit,         "phi", "link", "mean")
  test_dataframe_output(fit,         "phi", "link", "samples")
  test_dataframe_output(fit,         "phi", "response", "quantiles")
  test_dataframe_output(fit,         "phi", "response", "mean")
  test_dataframe_output(fit,         "phi", "response", "samples")
  test_dataframe_output(fit_unnamed, "phi", "link", "quantiles")
  test_dataframe_output(fit_unnamed, "phi", "link", "mean")
  test_dataframe_output(fit_unnamed, "phi", "link", "samples")
  test_dataframe_output(fit_unnamed, "phi", "response", "quantiles")
  test_dataframe_output(fit_unnamed, "phi", "response", "mean")
  test_dataframe_output(fit_unnamed, "phi", "response", "samples")

  # Tests for theta
  test_dataframe_output(fit,         "theta", "link", "quantiles")
  test_dataframe_output(fit,         "theta", "link", "mean")
  test_dataframe_output(fit,         "theta", "link", "samples")
  test_dataframe_output(fit,         "theta", "response", "quantiles")
  test_dataframe_output(fit,         "theta", "response", "mean")
  test_dataframe_output(fit,         "theta", "response", "samples")
  test_dataframe_output(fit_unnamed, "theta", "link", "quantiles")
  test_dataframe_output(fit_unnamed, "theta", "link", "mean")
  test_dataframe_output(fit_unnamed, "theta", "link", "samples")
  test_dataframe_output(fit_unnamed, "theta", "response", "quantiles")
  test_dataframe_output(fit_unnamed, "theta", "response", "mean")
  test_dataframe_output(fit_unnamed, "theta", "response", "samples")
})
