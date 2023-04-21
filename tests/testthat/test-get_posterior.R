### Tests for get_post_samples ------------------------------------------------
test_that("Extracted samples are correct when proper parameter names are given", {
    # TODO: test models with formula_* / formula_*_shared formulae.

    lpar <- c("z", "pi", "phi", "theta", "psi",
              "alpha", "beta", "gamma",
#              "alpha_share", "beta_shared", "gamma_shared",
              "Mu", "sigma", "rho")        

    I <- 2; J <- 2; K <- 2
    y <- array(sample.int(I * J * K), dim = c(I, J, K))
    spec_cov <- list(cov1 = rnorm(I))
    site_cov <- list(cov2 = rnorm(J), cov3 = factor(1:J))
    repl_cov <- list(cov4 = matrix(rnorm(J * K), J, K))
    data <- occumbData(
        y = y,
        spec_cov = spec_cov,
        site_cov = site_cov,
        repl_cov = repl_cov)

    fit <- occumb(data = data,
                  n.chains = 1, n.burnin = 10, n.thin = 1, n.iter = 20,
                  verbose = FALSE)

    for (i in seq_along(lpar)) {
        test <- get_post_samples(fit, lpar[i])
        attr(test, "dimension") <- NULL
        attr(test, "label") <- NULL
        expect_identical(test,
                         eval(parse(text = paste0("fit@fit$sims.list$", lpar[i]))))
    }
})

# TODO: Add tests for additional attributes supplied by get_post_samples().

