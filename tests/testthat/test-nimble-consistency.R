skip_on_cran()
skip_if_not_installed("nimble")
skip_if_not(identical(Sys.getenv("OCCUMB_TEST_NIMBLE_CONSISTENCY"), "true"),
            message = "Engine consistency tests skipped. Set OCCUMB_TEST_NIMBLE_CONSISTENCY=true to run.")

### Tests for engine consistency between NIMBLE and JAGS ----------------------
### These tests verify that occumbFit objects produced by the NIMBLE engine
### have a structure consistent with JAGS engine results, ensuring all
### downstream functions work equally with both engines.

# Null model fixtures
set.seed(42)
I <- 2; J <- 2; K <- 2
y <- array(sample.int(I * J * K), dim = c(I, J, K))
dimnames(y) <- list(c("sp1", "sp2"), c("s1", "s2"), c("r1", "r2"))
data_null <- occumbData(y = y)

mcmc_args <- list(n.chains = 2, n.burnin = 10, n.thin = 1, n.iter = 20,
                  verbose = FALSE)

fit_jags <- do.call(occumb, c(list(data = data_null), mcmc_args))
fit_nimble <- do.call(occumb,
                      c(list(data = data_null, engine = "NIMBLE"), mcmc_args))

# Covariate model fixtures
data_cov <- occumbData(
  y = y,
  spec_cov = list(cov1 = rnorm(I)),
  site_cov = list(cov2 = rnorm(J))
)

fit_jags_cov <- do.call(occumb,
  c(list(formula_psi_shared = ~ cov1, formula_theta = ~ cov2,
         data = data_cov), mcmc_args))
fit_nimble_cov <- do.call(occumb,
  c(list(formula_psi_shared = ~ cov1, formula_theta = ~ cov2,
         data = data_cov, engine = "NIMBLE"), mcmc_args))

null_model_params <- c("z", "pi", "phi", "theta", "psi",
                       "alpha", "beta", "gamma", "Mu", "sigma", "rho")

## Group 1: occumbFit slot-level structure ------------------------------------
test_that("occumbFit slots have correct types for both engines", {
  expect_identical(fit_jags@engine, "jags")
  expect_identical(fit_nimble@engine, "nimble")

  expect_s4_class(fit_jags@data, "occumbData")
  expect_s4_class(fit_nimble@data, "occumbData")

  expect_identical(names(fit_jags@occumb_args), names(fit_nimble@occumb_args))

  expect_true(inherits(fit_jags@fit, "jagsUI"))
  expect_true(inherits(fit_nimble@fit, "jagsUI"))
})

## Group 2: fit@fit top-level elements ----------------------------------------
test_that("fit@fit has required top-level elements for both engines", {
  required_names <- c("sims.list", "summary", "mcmc.info",
                      "parameters", "parallel", "run.date", "samples")
  expect_true(all(required_names %in% names(fit_jags@fit)))
  expect_true(all(required_names %in% names(fit_nimble@fit)))
})

## Group 3: $sims.list structure ----------------------------------------------
# Parameter order in sims.list may differ between engines; tests compare
# by name (set equality and name-based access) so order does not matter.
test_that("sims.list has the same parameters with consistent dimensions", {
  names_jags   <- names(fit_jags@fit$sims.list)
  names_nimble <- names(fit_nimble@fit$sims.list)

  # NIMBLE does not include "deviance" in sims.list
  names_jags_no_deviance <- setdiff(names_jags, "deviance")

  # Same set of parameter names (order-independent, excluding deviance)
  expect_setequal(names_jags_no_deviance, names_nimble)

  # Each parameter has the same array shape (excluding the samples dimension)
  for (param in names_jags_no_deviance) {
    dim_jags   <- dim(fit_jags@fit$sims.list[[param]])
    dim_nimble <- dim(fit_nimble@fit$sims.list[[param]])
    if (!is.null(dim_jags) && length(dim_jags) > 1) {
      expect_identical(dim_jags[-1], dim_nimble[-1],
                       info = paste("sims.list dim mismatch for", param))
    }
    # Both should be numeric
    expect_true(is.numeric(fit_jags@fit$sims.list[[param]]),
                info = paste("JAGS sims.list not numeric for", param))
    expect_true(is.numeric(fit_nimble@fit$sims.list[[param]]),
                info = paste("NIMBLE sims.list not numeric for", param))
  }
})

## Group 4: $summary structure ------------------------------------------------
test_that("summary matrix has consistent structure", {
  expected_cols <- c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%",
                     "Rhat", "n.eff", "overlap0", "f")
  expect_identical(colnames(fit_jags@fit$summary), expected_cols)
  expect_identical(colnames(fit_nimble@fit$summary), expected_cols)

  # NIMBLE does not include "deviance" in summary
  rownames_jags_no_deviance <- setdiff(rownames(fit_jags@fit$summary), "deviance")
  expect_setequal(rownames_jags_no_deviance,
                  rownames(fit_nimble@fit$summary))
})

## Group 5: $mcmc.info structure ----------------------------------------------
test_that("mcmc.info has required fields for NIMBLE", {
  required_fields <- c("n.chains", "n.iter", "n.burnin", "n.thin",
                       "n.samples", "elapsed.mins")
  expect_true(all(required_fields %in% names(fit_nimble@fit$mcmc.info)))

  expect_true(is.numeric(fit_nimble@fit$mcmc.info$n.chains))
  expect_true(is.numeric(fit_nimble@fit$mcmc.info$n.samples))
  expect_true(is.numeric(fit_nimble@fit$mcmc.info$elapsed.mins))
})

## Group 6: $parameters, $parallel, $run.date, $samples -----------------------
test_that("Auxiliary fit elements are consistent", {
  expect_true(is.character(fit_nimble@fit$parameters))
  # NIMBLE does not include "deviance" in parameters
  expect_setequal(setdiff(fit_jags@fit$parameters, "deviance"),
                  fit_nimble@fit$parameters)

  expect_true(is.logical(fit_nimble@fit$parallel))

  expect_true(inherits(fit_nimble@fit$run.date, "POSIXt"))

  expect_true(inherits(fit_nimble@fit$samples, "mcmc.list"))
  expect_equal(coda::nchain(fit_nimble@fit$samples),
               fit_nimble@fit$mcmc.info$n.chains)
})

## Group 7: Downstream functions work without error ---------------------------
test_that("show() works for both engines", {
  expect_output(show(fit_jags))
  expect_output(show(fit_nimble))
})

test_that("summary() works for both engines", {
  expect_output(summary(fit_jags))
  expect_output(summary(fit_nimble))
})

test_that("plot() works for both engines", {
  pdf(nullfile())
  on.exit(dev.off())
  expect_no_error(plot(fit_jags))
  expect_no_error(plot(fit_nimble))
})

test_that("get_post_samples() works for both engines on all parameters", {
  for (param in null_model_params) {
    samples_jags   <- get_post_samples(fit_jags, param)
    samples_nimble <- get_post_samples(fit_nimble, param)
    if (!is.null(dim(samples_jags)) && length(dim(samples_jags)) > 1) {
      expect_identical(dim(samples_jags)[-1], dim(samples_nimble)[-1],
                       info = paste("get_post_samples dim mismatch for", param))
    }
  }
})

test_that("get_post_summary() works for both engines on all parameters", {
  for (param in null_model_params) {
    summary_jags   <- get_post_summary(fit_jags, param)
    summary_nimble <- get_post_summary(fit_nimble, param)
    if (!is.null(dim(summary_jags))) {
      expect_identical(colnames(summary_jags), colnames(summary_nimble),
                       info = paste("get_post_summary colnames mismatch for", param))
    } else {
      expect_identical(names(summary_jags), names(summary_nimble),
                       info = paste("get_post_summary names mismatch for", param))
    }
  }
})

test_that("predict() works for both engines", {
  for (param in c("phi", "theta", "psi")) {
    expect_no_error(predict(fit_jags, parameter = param))
    expect_no_error(predict(fit_nimble, parameter = param))
  }
})

test_that("gof() works for both engines", {
  expect_no_error(gof(fit_jags, plot = FALSE))
  expect_no_error(gof(fit_nimble, plot = FALSE))
})

## Group 8: Covariate model consistency ---------------------------------------
test_that("Covariate model produces consistent structure", {
  expect_true("gamma_shared" %in% names(fit_jags_cov@fit$sims.list))
  expect_true("gamma_shared" %in% names(fit_nimble_cov@fit$sims.list))

  # NIMBLE does not include "deviance" in sims.list or summary
  expect_setequal(setdiff(names(fit_jags_cov@fit$sims.list), "deviance"),
                  names(fit_nimble_cov@fit$sims.list))
  expect_setequal(setdiff(rownames(fit_jags_cov@fit$summary), "deviance"),
                  rownames(fit_nimble_cov@fit$summary))
})

test_that("Downstream functions work on covariate model for both engines", {
  expect_no_error(get_post_samples(fit_jags_cov, "gamma_shared"))
  expect_no_error(get_post_samples(fit_nimble_cov, "gamma_shared"))

  for (param in c("phi", "theta", "psi")) {
    expect_no_error(predict(fit_jags_cov, parameter = param))
    expect_no_error(predict(fit_nimble_cov, parameter = param))
  }
})
