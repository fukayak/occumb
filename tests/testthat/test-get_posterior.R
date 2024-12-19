### Test data ------------------------------------------------------------------
I <- 2
J <- 3
K <- 4
cov1 <- rnorm(I)
cov2 <- factor(1:I)
cov3 <- rnorm(J)
cov4 <- factor(1:J)
cov5 <- matrix(rnorm(J * K), nrow = J)
cov6 <- matrix(rep(factor(1:K), each = J), nrow = J)
y <- array(sample.int(1E3, I * J * K), dim = c(I, J, K))
y[1, 1, ] <- 0
y[1, 2, ] <- 0
data <- occumbData(y = y,
                   spec_cov = list(cov1 = cov1, cov2 = cov2),
                   site_cov = list(cov3 = cov3, cov4 = cov4),
                   repl_cov = list(cov5 = cov5, cov6 = cov6))

y_named <- array(sample.int(1E3, I * J * K), dim = c(I, J, K))
y_named[1, 1, ] <- 0
y_named[1, 2, ] <- 0
dimnames(y_named) <- list(sprintf("sp. %s", 1:2),
                          sprintf("site %s", 1:3),
                          sprintf("repl. %s", 1:4))
data_named <- occumbData(y = y_named,
                         spec_cov = list(cov1 = cov1, cov2 = cov2),
                         site_cov = list(cov3 = cov3, cov4 = cov4),
                         repl_cov = list(cov5 = cov5, cov6 = cov6))

fit1 <- occumb(formula_phi = ~ cov6,
               formula_theta = ~ cov4,
               formula_psi = ~ 1,
               formula_phi_shared = ~ 1,
               formula_theta_shared = ~ cov2 * cov3,
               formula_psi_shared = ~ cov2,
               data = data,
               n.chains = 1, n.burnin = 5, n.thin = 1, n.iter = 10,
               verbose = FALSE)

fit2 <- occumb(formula_phi = ~ cov6,
               formula_theta = ~ cov4,
               formula_psi = ~ 1,
               formula_phi_shared = ~ 1,
               formula_theta_shared = ~ cov2 * cov3,
               formula_psi_shared = ~ cov2,
               data = data_named,
               n.chains = 1, n.burnin = 5, n.thin = 1, n.iter = 10,
               verbose = FALSE)

fit3 <- occumb(formula_phi = ~ 1,
               formula_theta = ~ 1,
               formula_psi = ~ 1,
               formula_phi_shared = ~ cov6,
               formula_theta_shared = ~ cov4,
               formula_psi_shared = ~ cov3,
               data = data_named,
               n.chains = 1, n.burnin = 5, n.thin = 1, n.iter = 10,
               verbose = FALSE)

lpar <- c("z", "pi", "phi", "theta", "psi",
          "alpha", "beta", "gamma",
          "alpha_shared", "beta_shared", "gamma_shared",
          "Mu", "sigma", "rho")

### Tests for get_post_samples ------------------------------------------------
test_that("Extracted samples and attributes are correct when proper parameter names are given", {

  ## Tests for z
  i <- 1
  # Unnamed y
  test <- get_post_samples(fit1, lpar[i])
  expect_identical(attributes(test)$dimension, c("Sample", "Species", "Site"))
  attr(test, "dimension") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit1@fit$sims.list$", lpar[i]))))
  # Named y
  test <- get_post_samples(fit2, lpar[i])
  expect_identical(attributes(test)$dimension, c("Sample", "Species", "Site"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Species = dimnames(data_named@y)[[1]],
                        Site = dimnames(data_named@y)[[2]]))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit2@fit$sims.list$", lpar[i]))))

  ## Tests for pi
  i <- 2
  # Unnamed y
  test <- get_post_samples(fit1, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Site", "Replicate"))
  attr(test, "dimension") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit1@fit$sims.list$", lpar[i]))))
  # Named y
  test <- get_post_samples(fit2, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Site", "Replicate"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Species = dimnames(data_named@y)[[1]],
                        Site = dimnames(data_named@y)[[2]],
                        Replicate = dimnames(data_named@y)[[3]]))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit2@fit$sims.list$", lpar[i]))))

  ## Tests for phi (= "ijk")
  i <- 3
  # Unnamed y
  test <- get_post_samples(fit1, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Site", "Replicate"))
  attr(test, "dimension") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit1@fit$sims.list$", lpar[i]))))
  # Named y
  test <- get_post_samples(fit2, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Site", "Replicate"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Species = dimnames(data_named@y)[[1]],
                        Site = dimnames(data_named@y)[[2]],
                        Replicate = dimnames(data_named@y)[[3]]))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit2@fit$sims.list$", lpar[i]))))

  ## Tests for theta (= "ij")
  i <- 4
  # Unnamed y
  test <- get_post_samples(fit1, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Site"))
  attr(test, "dimension") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit1@fit$sims.list$", lpar[i]))))
  # Named y
  test <- get_post_samples(fit2, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Site"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Species = dimnames(data_named@y)[[1]],
                        Site = dimnames(data_named@y)[[2]]))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit2@fit$sims.list$", lpar[i]))))

  ## Tests for psi (= "i")
  i <- 5
  # Unnamed y
  test <- get_post_samples(fit1, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species"))
  attr(test, "dimension") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit1@fit$sims.list$", lpar[i]))))
  # Named y
  test <- get_post_samples(fit2, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Species = dimnames(data_named@y)[[1]]))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit2@fit$sims.list$", lpar[i]))))

  ## Tests for alpha (~ cov6)
  i <- 6
  # Unnamed y
  test <- get_post_samples(fit1, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Species = NULL,
                        Effects = c("(Intercept)", "cov62", "cov63", "cov64")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit1@fit$sims.list$", lpar[i]))))
  # Named y
  test <- get_post_samples(fit2, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Species = dimnames(data_named@y)[[1]],
                        Effects = c("(Intercept)", "cov62", "cov63", "cov64")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit2@fit$sims.list$", lpar[i]))))

  ## Tests for beta (~ cov4)
  i <- 7
  # Unnamed y
  test <- get_post_samples(fit1, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Species = NULL,
                        Effects = c("(Intercept)", "cov42", "cov43")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit1@fit$sims.list$", lpar[i]))))
  # Named y
  test <- get_post_samples(fit2, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Species = dimnames(data_named@y)[[1]],
                        Effects = c("(Intercept)", "cov42", "cov43")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit2@fit$sims.list$", lpar[i]))))

  ## Tests for gamma (~ 1)
  i <- 8
  # Unnamed y
  test <- get_post_samples(fit1, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Species = NULL,
                        Effects = "(Intercept)"))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit1@fit$sims.list$", lpar[i]))))
  # Named y
  test <- get_post_samples(fit2, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Species = dimnames(data_named@y)[[1]],
                        Effects = "(Intercept)"))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit2@fit$sims.list$", lpar[i]))))

  ## Tests for alpha_shared (~ 1)
  i <- 9
  # Unnamed y
  expect_error(get_post_samples(fit1, lpar[i]),
               "alpha_shared is not included in the fitted model")
  # Named y
  expect_error(get_post_samples(fit2, lpar[i]),
               "alpha_shared is not included in the fitted model")

  ## Tests for beta_shared (~ cov2 * cov3)
  i <- 10
  # Unnamed y
  test <- get_post_samples(fit1, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Effects = c("cov22", "cov3", "cov22:cov3")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit1@fit$sims.list$", lpar[i]))))
  # Named y
  test <- get_post_samples(fit2, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Effects = c("cov22", "cov3", "cov22:cov3")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit2@fit$sims.list$", lpar[i]))))

  ## Tests for gamma_shared (~ cov2)
  i <- 11
  # Unnamed y
  test <- get_post_samples(fit1, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Effects = c("cov22")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit1@fit$sims.list$", lpar[i]))))
  # Named y
  test <- get_post_samples(fit2, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Effects = c("cov22")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit2@fit$sims.list$", lpar[i]))))

  ## Tests for Mu/sigma
  for (i in 12:13) {
    # Unnamed y
    test <- get_post_samples(fit1, lpar[i])
    expect_identical(attributes(test)$dimension,
                     c("Sample", "Effects"))
    expect_identical(attributes(test)$label,
                     list(Sample = NULL,
                          Effects = c("phi | (Intercept)", "phi | cov62",
                                      "phi | cov63", "phi | cov64",
                                      "theta | (Intercept)", "theta | cov42", "theta | cov43",
                                      "psi | (Intercept)")))
    attr(test, "dimension") <- NULL
    attr(test, "label") <- NULL
    expect_identical(test,
                     eval(parse(text = paste0("fit1@fit$sims.list$", lpar[i]))))
    # Named y
    test <- get_post_samples(fit2, lpar[i])
    expect_identical(attributes(test)$dimension,
                     c("Sample", "Effects"))
    expect_identical(attributes(test)$label,
                     list(Sample = NULL,
                          Effects = c("phi | (Intercept)", "phi | cov62",
                                      "phi | cov63", "phi | cov64",
                                      "theta | (Intercept)", "theta | cov42", "theta | cov43",
                                      "psi | (Intercept)")))
    attr(test, "dimension") <- NULL
    attr(test, "label") <- NULL
    expect_identical(test,
                     eval(parse(text = paste0("fit2@fit$sims.list$", lpar[i]))))
  }

  ## Tests for rho
  i <- 14
  # Unnamed y
  test <- get_post_samples(fit1, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Effects 1", "Effects 2"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Effects1 = c("phi | (Intercept)", "phi | cov62",
                                     "phi | cov63", "phi | cov64",
                                     "theta | (Intercept)", "theta | cov42", "theta | cov43"),
                        Effects2 = c("phi | (Intercept)", "phi | cov62",
                                     "phi | cov63", "phi | cov64",
                                     "theta | (Intercept)", "theta | cov42", "theta | cov43",
                                     "psi | (Intercept)")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- eval(parse(text = paste0("fit1@fit$sims.list$", lpar[i])))
  expect_identical(test, ans)
  # Named y
  test <- get_post_samples(fit2, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Effects 1", "Effects 2"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Effects1 = c("phi | (Intercept)", "phi | cov62",
                                     "phi | cov63", "phi | cov64",
                                     "theta | (Intercept)", "theta | cov42", "theta | cov43"),
                        Effects2 = c("phi | (Intercept)", "phi | cov62",
                                     "phi | cov63", "phi | cov64",
                                     "theta | (Intercept)", "theta | cov42", "theta | cov43",
                                     "psi | (Intercept)")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- eval(parse(text = paste0("fit2@fit$sims.list$", lpar[i])))
  expect_identical(test, ans)
})

test_that("$label$Effects attributes for alpha/beta/gamma are given correctly when site_cov or repl_cov is specified for shared parameters", {
  ## Tests for alpha
  i <- 6
  test <- get_post_samples(fit3, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Species = dimnames(data_named@y)[[1]],
                        Effects = "(Intercept)"))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit3@fit$sims.list$", lpar[i]))))

  ## Tests for beta
  i <- 7
  test <- get_post_samples(fit3, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Species = dimnames(data_named@y)[[1]],
                        Effects = "(Intercept)"))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit3@fit$sims.list$", lpar[i]))))

  ## Tests for gamma
  i <- 8
  test <- get_post_samples(fit3, lpar[i])
  expect_identical(attributes(test)$dimension,
                   c("Sample", "Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Sample = NULL,
                        Species = dimnames(data_named@y)[[1]],
                        Effects = "(Intercept)"))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  expect_identical(test,
                   eval(parse(text = paste0("fit3@fit$sims.list$", lpar[i]))))
})

### Tests for get_post_summary -------------------------------------------------
test_that("Extracted tables and attributes are correct when proper parameter names are given", {

  ## Tests for z
  i <- 1
  # Unnamed y
  test <- get_post_summary(fit1, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Site"))
  attr(test, "dimension") <- NULL
  ans <- fit1@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit1@fit$summary)), ]
  expect_identical(test, ans)
  # Named y
  test <- get_post_summary(fit2, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Site"))
  expect_identical(attributes(test)$label,
                   list(Species = dimnames(data_named@y)[[1]],
                        Site = dimnames(data_named@y)[[2]]))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit2@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit2@fit$summary)), ]
  expect_identical(test, ans)

  ## Tests for pi
  i <- 2
  # Unnamed y
  test <- get_post_summary(fit1, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Site", "Replicate"))
  attr(test, "dimension") <- NULL
  ans <- fit1@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit1@fit$summary)), ]
  expect_identical(test, ans)
  # Named y
  test <- get_post_summary(fit2, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Site", "Replicate"))
  expect_identical(attributes(test)$label,
                   list(Species = dimnames(data_named@y)[[1]],
                        Site = dimnames(data_named@y)[[2]],
                        Replicate = dimnames(data_named@y)[[3]]))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit2@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit2@fit$summary)), ]
  expect_identical(test, ans)

  ## Tests for phi (= "ijk")
  i <- 3
  # Unnamed y
  test <- get_post_summary(fit1, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Site", "Replicate"))
  attr(test, "dimension") <- NULL
  ans <- fit1@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit1@fit$summary)), ]
  expect_identical(test, ans)
  # Named y
  test <- get_post_summary(fit2, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Site", "Replicate"))
  expect_identical(attributes(test)$label,
                   list(Species = dimnames(data_named@y)[[1]],
                        Site = dimnames(data_named@y)[[2]],
                        Replicate = dimnames(data_named@y)[[3]]))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit2@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit2@fit$summary)), ]
  expect_identical(test, ans)

  ## Tests for theta (= "ij")
  i <- 4
  # Unnamed y
  test <- get_post_summary(fit1, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Site"))
  attr(test, "dimension") <- NULL
  ans <- fit1@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit1@fit$summary)), ]
  expect_identical(test, ans)
  # Named y
  test <- get_post_summary(fit2, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Site"))
  expect_identical(attributes(test)$label,
                   list(Species = dimnames(data_named@y)[[1]],
                        Site = dimnames(data_named@y)[[2]]))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit2@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit2@fit$summary)), ]
  expect_identical(test, ans)

  ## Tests for psi (= "i")
  i <- 5
  # Unnamed y
  test <- get_post_summary(fit1, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species"))
  attr(test, "dimension") <- NULL
  ans <- fit1@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit1@fit$summary)), ]
  expect_identical(test, ans)
  # Named y
  test <- get_post_summary(fit2, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species"))
  expect_identical(attributes(test)$label,
                   list(Species = dimnames(data_named@y)[[1]]))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit2@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit2@fit$summary)), ]
  expect_identical(test, ans)

  ## Tests for alpha (~ cov6)
  i <- 6
  # Unnamed y
  test <- get_post_summary(fit1, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Species = NULL,
                        Effects = c("(Intercept)", "cov62", "cov63", "cov64")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit1@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit1@fit$summary)), ]
  expect_identical(test, ans)
  # Named y
  test <- get_post_summary(fit2, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Species = dimnames(data_named@y)[[1]],
                        Effects = c("(Intercept)", "cov62", "cov63", "cov64")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit2@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit2@fit$summary)), ]
  expect_identical(test, ans)

  ## Tests for beta (~ cov4)
  i <- 7
  # Unnamed y
  test <- get_post_summary(fit1, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Species = NULL,
                        Effects = c("(Intercept)", "cov42", "cov43")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit1@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit1@fit$summary)), ]
  expect_identical(test, ans)
  # Named y
  test <- get_post_summary(fit2, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Species = dimnames(data_named@y)[[1]],
                        Effects = c("(Intercept)", "cov42", "cov43")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit2@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit2@fit$summary)), ]
  expect_identical(test, ans)

  ## Tests for gamma (~ 1)
  i <- 8
  # Unnamed y
  test <- get_post_summary(fit1, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Species = NULL,
                        Effects = "(Intercept)"))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit1@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit1@fit$summary)), ]
  expect_identical(test, ans)
  # Named y
  test <- get_post_summary(fit2, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Species = dimnames(data_named@y)[[1]],
                        Effects = "(Intercept)"))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit2@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit2@fit$summary)), ]
  expect_identical(test, ans)

  ## Tests for alpha_shared (~ 1)
  i <- 9
  # Unnamed y
  expect_error(get_post_samples(fit1, lpar[i]),
               "alpha_shared is not included in the fitted model")
  # Named y
  expect_error(get_post_samples(fit2, lpar[i]),
               "alpha_shared is not included in the fitted model")

  ## Tests for beta_shared (~ cov2 * cov3)
  i <- 10
  # Unnamed y
  test <- get_post_summary(fit1, lpar[i])
  expect_identical(attributes(test)$dimension, c("Effects"))
  expect_identical(attributes(test)$label,
                   list(Effects = c("cov22", "cov3", "cov22:cov3")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit1@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit1@fit$summary)), ]
  expect_identical(test, ans)
  # Named y
  test <- get_post_summary(fit2, lpar[i])
  expect_identical(attributes(test)$dimension, c("Effects"))
  expect_identical(attributes(test)$label,
                   list(Effects = c("cov22", "cov3", "cov22:cov3")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit2@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit2@fit$summary)), ]
  expect_identical(test, ans)

  ## Tests for gamma_shared (~ cov2)
  i <- 11
  # Unnamed y
  test <- get_post_summary(fit1, lpar[i])
  expect_identical(attributes(test)$dimension, c("Effects"))
  expect_identical(attributes(test)$label,
                   list(Effects = "cov22"))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit1@fit$summary["gamma_shared", ]
  expect_identical(test, ans)
  # Named y
  test <- get_post_summary(fit2, lpar[i])
  expect_identical(attributes(test)$dimension, c("Effects"))
  expect_identical(attributes(test)$label,
                   list(Effects = "cov22"))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit2@fit$summary["gamma_shared", ]
  expect_identical(test, ans)

  ## Tests for Mu/sigma
  for (i in 12:13) {
    # Unnamed y
    test <- get_post_summary(fit1, lpar[i])
    expect_identical(attributes(test)$dimension, c("Effects"))
    expect_identical(attributes(test)$label,
                     list(Effects = c("phi | (Intercept)", "phi | cov62",
                                      "phi | cov63", "phi | cov64",
                                      "theta | (Intercept)", "theta | cov42", "theta | cov43",
                                      "psi | (Intercept)")))
    attr(test, "dimension") <- NULL
    attr(test, "label") <- NULL
    ans <- fit1@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit1@fit$summary)), ]
    expect_identical(test, ans)
    # Named y
    test <- get_post_summary(fit2, lpar[i])
    expect_identical(attributes(test)$dimension, c("Effects"))
    expect_identical(attributes(test)$label,
                     list(Effects = c("phi | (Intercept)", "phi | cov62",
                                      "phi | cov63", "phi | cov64",
                                      "theta | (Intercept)", "theta | cov42", "theta | cov43",
                                      "psi | (Intercept)")))
    attr(test, "dimension") <- NULL
    attr(test, "label") <- NULL
    ans <- fit2@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit2@fit$summary)), ]
    expect_identical(test, ans)
  }

  ## Tests for rho
  i <- 14
  # Unnamed y
  test <- get_post_summary(fit1, lpar[i])
  expect_identical(attributes(test)$dimension, c("Effects 1", "Effects 2"))
  expect_identical(attributes(test)$label,
                   list(Effects1 = c("phi | (Intercept)", "phi | cov62",
                                     "phi | cov63", "phi | cov64",
                                     "theta | (Intercept)", "theta | cov42", "theta | cov43"),
                        Effects2 = c("phi | (Intercept)", "phi | cov62",
                                     "phi | cov63", "phi | cov64",
                                     "theta | (Intercept)", "theta | cov42", "theta | cov43",
                                     "psi | (Intercept)")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit1@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit1@fit$summary)), ]
  expect_identical(test, ans)
  # Named y
  test <- get_post_summary(fit2, lpar[i])
  expect_identical(attributes(test)$dimension, c("Effects 1", "Effects 2"))
  expect_identical(attributes(test)$label,
                   list(Effects1 = c("phi | (Intercept)", "phi | cov62",
                                     "phi | cov63", "phi | cov64",
                                     "theta | (Intercept)", "theta | cov42", "theta | cov43"),
                        Effects2 = c("phi | (Intercept)", "phi | cov62",
                                     "phi | cov63", "phi | cov64",
                                     "theta | (Intercept)", "theta | cov42", "theta | cov43",
                                     "psi | (Intercept)")))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit2@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit2@fit$summary)), ]
  expect_identical(test, ans)

})
test_that("$label$Effects attributes for alpha/beta/gamma are given correctly when site_cov or repl_cov is specified for shared parameters", {
  ## Tests for alpha
  i <- 6
  test <- get_post_summary(fit3, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Species = dimnames(data_named@y)[[1]],
                        Effects = "(Intercept)"))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit3@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit3@fit$summary)), ]
  expect_identical(test, ans)

  ## Tests for beta
  i <- 7
  test <- get_post_summary(fit3, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Species = dimnames(data_named@y)[[1]],
                        Effects = "(Intercept)"))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit3@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit3@fit$summary)), ]
  expect_identical(test, ans)

  ## Tests for gamma
  i <- 8
  test <- get_post_summary(fit3, lpar[i])
  expect_identical(attributes(test)$dimension, c("Species", "Effects"))
  expect_identical(attributes(test)$label,
                   list(Species = dimnames(data_named@y)[[1]],
                        Effects = "(Intercept)"))
  attr(test, "dimension") <- NULL
  attr(test, "label") <- NULL
  ans <- fit3@fit$summary[grep(paste0(lpar[i], "\\["), rownames(fit3@fit$summary)), ]
  expect_identical(test, ans)
})
