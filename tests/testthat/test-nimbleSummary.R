### Unit tests for nimbleSummary.R ---------------------------------------------
### These tests verify the internal functions that convert NIMBLE MCMC output
### into jagsUI-compatible summary objects, using hand-crafted mcmc.list
### fixtures (no NIMBLE/JAGS dependency required).

# --- Test fixtures -----------------------------------------------------------
make_test_mcmc <- function(n_iter = 20, n_chain = 2, seed = 42) {
  set.seed(seed)
  params <- c("alpha", "beta[1]", "beta[2]", "gamma[1,1]", "gamma[1,2]")
  chains <- lapply(seq_len(n_chain), function(ch) {
    mat <- matrix(rnorm(n_iter * length(params)), nrow = n_iter)
    colnames(mat) <- params
    coda::as.mcmc(mat)
  })
  coda::as.mcmc.list(chains)
}

make_single_param_mcmc <- function(name = "mu", n_iter = 20, n_chain = 2,
                                   seed = 42, fill = NULL) {
  set.seed(seed)
  chains <- lapply(seq_len(n_chain), function(ch) {
    if (is.null(fill)) {
      mat <- matrix(rnorm(n_iter), nrow = n_iter)
    } else {
      mat <- matrix(fill, nrow = n_iter, ncol = 1)
    }
    colnames(mat) <- name
    coda::as.mcmc(mat)
  })
  coda::as.mcmc.list(chains)
}

add_print_fields <- function(result) {
  result$modfile <- "test_model"
  result$mcmc.info$n.burnin <- 0
  result$mcmc.info$n.thin <- 1
  result$mcmc.info$elapsed.mins <- 0.1
  result$parallel <- FALSE
  result$run.date <- Sys.time()
  result
}

mcmc2 <- make_test_mcmc()
mcmc1 <- make_test_mcmc(n_chain = 1)
single_param_mcmc <- make_single_param_mcmc()
stats2 <- occumb:::calc_stats(mcmc2)
stats_coda <- occumb:::calc_stats(mcmc2, coda_only = "alpha")
summary2 <- occumb:::nimbleSummary(mcmc2)
summary1 <- occumb:::nimbleSummary(mcmc1)

## Low-level helper functions --------------------------------------------------

test_that("strip_params removes brackets", {
  expect_equal(
    occumb:::strip_params(c("alpha", "beta[1]", "beta[2]", "gamma[1,1]")),
    c("alpha", "beta", "beta", "gamma")
  )
})

test_that("strip_params with unique = TRUE returns unique base names", {
  expect_equal(
    occumb:::strip_params(c("alpha", "beta[1]", "beta[2]"), unique = TRUE),
    c("alpha", "beta")
  )
})

test_that("which_params returns correct indices for array parameters", {
  raw <- c("alpha", "beta[1]", "beta[2]", "gamma[1,1]")
  expect_equal(occumb:::which_params("beta", raw), c(2L, 3L))
  expect_equal(occumb:::which_params("alpha", raw), 1L)
  expect_null(occumb:::which_params("nonexistent", raw))
})

test_that("param_names extracts variable names from mcmc.list", {
  expect_equal(
    occumb:::param_names(mcmc2),
    c("alpha", "beta[1]", "beta[2]", "gamma[1,1]", "gamma[1,2]")
  )
  expect_equal(
    occumb:::param_names(mcmc2, simplify = TRUE),
    c("alpha", "beta", "gamma")
  )
})

test_that("match_params matches scalar and array parameters", {
  raw <- c("alpha", "beta[1]", "beta[2]", "gamma[1,1]")
  expect_equal(occumb:::match_params("alpha", raw), "alpha")
  expect_equal(occumb:::match_params("beta", raw), c("beta[1]", "beta[2]"))
  expect_equal(occumb:::match_params("gamma", raw), "gamma[1,1]")
  expect_length(occumb:::match_params("nonexistent", raw), 0)
  expect_equal(
    occumb:::match_params(c("beta", "alpha"), raw),
    c("beta[1]", "beta[2]", "alpha")
  )
})

test_that("get_inds extracts index matrix from parameter names", {
  inds <- occumb:::get_inds("beta", c("beta[1]", "beta[2]", "beta[3]"))
  expect_equal(inds, matrix(1:3, ncol = 1))
  inds2 <- occumb:::get_inds("gamma", c("gamma[1,1]", "gamma[1,2]", "gamma[2,1]"))
  expect_equal(inds2, matrix(c(1, 1, 1, 2, 2, 1), ncol = 2, byrow = TRUE))
})

test_that("fill_array fills correct positions and leaves others NA", {
  indices <- matrix(c(1, 2, 2, 1), ncol = 2, byrow = TRUE)
  result <- occumb:::fill_array(c(10, 20), indices)
  expect_equal(dim(result), c(2, 2))
  expect_equal(result[1, 2], 10)
  expect_equal(result[2, 1], 20)
  expect_true(is.na(result[1, 1]))
  expect_true(is.na(result[2, 2]))
})

test_that("has_one_parameter returns TRUE for single-column mcmc.list", {
  expect_true(occumb:::has_one_parameter(single_param_mcmc))
  expect_false(occumb:::has_one_parameter(mcmc2))
})

## Statistical computation functions -------------------------------------------

test_that("overlap_0 detects whether 0 is in the interval", {
  expect_equal(occumb:::overlap_0(-1, 1), 1)
  expect_equal(occumb:::overlap_0(0.5, 2), 0)
  expect_equal(occumb:::overlap_0(-3, -0.1), 0)
})

test_that("calc_f returns proportion with same sign as mean", {
  expect_equal(occumb:::calc_f(c(1, 2, 3, 4), 2.5), 1)
  expect_equal(occumb:::calc_f(c(-1, -2, -3, -4), -2.5), 1)
  expect_equal(occumb:::calc_f(c(-1, 1, 2, 3), 1.25), 0.75)
})

test_that("calc_Rhat returns NA for single chain", {
  single <- make_single_param_mcmc(n_chain = 1)
  expect_true(is.na(occumb:::calc_Rhat(single)))
})

test_that("calc_Rhat returns a numeric value for multiple chains", {
  rhat <- occumb:::calc_Rhat(single_param_mcmc)
  expect_true(is.numeric(rhat))
  expect_true(is.finite(rhat))
})

test_that("calc_neff returns a positive integer for normal input", {
  multi <- make_single_param_mcmc(n_chain = 2, n_iter = 100)
  neff <- occumb:::calc_neff(multi)
  expect_true(is.numeric(neff))
  expect_true(neff >= 1)
})

test_that("calc_neff returns 1 for degenerate (zero-variance) input", {
  degenerate <- make_single_param_mcmc(name = "const", fill = 5)
  expect_equal(occumb:::calc_neff(degenerate), 1)
})

test_that("mcmc_to_mat converts mcmc.list to correct matrix shape", {
  single <- make_single_param_mcmc(n_iter = 10, n_chain = 3)
  mat <- occumb:::mcmc_to_mat(single)
  expect_equal(dim(mat), c(10, 3))
})

test_that("calc_param_stats returns 11 named statistics", {
  stats <- occumb:::calc_param_stats(single_param_mcmc, coda_only = FALSE)
  expected_names <- c("mean", "sd", "q2.5", "q25", "q50", "q75", "q97.5",
                      "overlap0", "f", "Rhat", "n.eff")
  expect_equal(names(stats), expected_names)
  expect_true(all(is.numeric(stats)))
})

test_that("calc_param_stats returns NA vector for all-Inf input", {
  inf_mcmc <- make_single_param_mcmc(name = "bad", fill = Inf)
  stats <- occumb:::calc_param_stats(inf_mcmc, coda_only = FALSE)
  expect_true(all(is.na(stats)))
})

test_that("calc_param_stats returns NA vector for all-NA input", {
  na_mcmc <- make_single_param_mcmc(name = "bad", fill = NA_real_)
  stats <- occumb:::calc_param_stats(na_mcmc, coda_only = FALSE)
  expect_true(all(is.na(stats)))
})

test_that("calc_param_stats with coda_only returns only mean", {
  stats <- occumb:::calc_param_stats(single_param_mcmc, coda_only = TRUE)
  expect_false(is.na(stats["mean"]))
  expect_true(all(is.na(stats[setdiff(names(stats), "mean")])))
})

## Aggregation and output formatting functions ---------------------------------

test_that("calc_stats produces correct dimensions for mixed parameters", {
  expect_equal(nrow(stats2), 5)
  expect_equal(ncol(stats2), 11)
  expect_equal(rownames(stats2),
               c("alpha", "beta[1]", "beta[2]", "gamma[1,1]", "gamma[1,2]"))
})

test_that("stat_summary_table renames quantile columns and reorders", {
  tbl <- occumb:::stat_summary_table(stats2, coda_only = NULL)
  expect_equal(
    colnames(tbl),
    c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%",
      "Rhat", "n.eff", "overlap0", "f")
  )
})

test_that("stat_summary_table excludes coda_only parameters", {
  tbl <- occumb:::stat_summary_table(stats_coda, coda_only = "alpha")
  expect_false("alpha" %in% rownames(tbl))
  expect_equal(nrow(tbl), 4)
})

test_that("get_posterior_array returns vector for scalar parameter", {
  result <- occumb:::get_posterior_array("alpha", mcmc2)
  expect_true(is.vector(result))
  expect_length(result, 40)
})

test_that("get_posterior_array returns correct array for array parameter", {
  result <- occumb:::get_posterior_array("beta", mcmc2)
  expect_equal(dim(result), c(40, 2))
})

test_that("get_posterior_array returns correct array for 2D parameter", {
  result <- occumb:::get_posterior_array("gamma", mcmc2)
  expect_equal(dim(result), c(40, 1, 2))
})

test_that("sims_list produces named list for all base parameters", {
  sl <- occumb:::sims_list(mcmc2)
  expect_equal(names(sl), c("alpha", "beta", "gamma"))
  expect_true(is.vector(sl$alpha))
  expect_equal(dim(sl$beta), c(40, 2))
})

test_that("get_stat_array returns scalar for scalar parameter", {
  result <- occumb:::get_stat_array("alpha", "mean", stats2)
  expect_length(result, 1)
  expect_true(is.numeric(result))
})

test_that("get_stat_array returns correct array for array parameter", {
  result <- occumb:::get_stat_array("beta", "mean", stats2)
  expect_length(result, 2)
})

test_that("all_stat_arrays converts overlap0 to logical", {
  arrays <- occumb:::all_stat_arrays(stats2, coda_only = NULL)
  expect_true(is.logical(arrays$overlap0$alpha))
  expect_true(is.logical(arrays$overlap0$beta))
})

test_that("all_stat_arrays has all expected stat names", {
  arrays <- occumb:::all_stat_arrays(stats2, coda_only = NULL)
  expected <- c("mean", "sd", "q2.5", "q25", "q50", "q75", "q97.5",
                "overlap0", "f", "Rhat", "n.eff")
  expect_equal(names(arrays), expected)
})

test_that("all_stat_arrays returns NA for non-mean stats of coda_only params", {
  arrays <- occumb:::all_stat_arrays(stats_coda, coda_only = "alpha")
  expect_true(is.numeric(arrays$mean$alpha))
  expect_true(is.na(arrays$sd$alpha))
  expect_true(is.na(arrays$Rhat$alpha))
  expect_true(is.na(arrays$n.eff$alpha))
})

test_that("calc_DIC returns NULL when DIC = FALSE", {
  expect_null(occumb:::calc_DIC(mcmc2, DIC = FALSE))
})

test_that("calc_DIC returns NULL when deviance is not in samples", {
  expect_null(occumb:::calc_DIC(mcmc2, DIC = TRUE))
})

test_that("calc_DIC computes pD and DIC when deviance is present", {
  set.seed(1)
  chains <- lapply(1:2, function(ch) {
    mat <- matrix(rexp(20, rate = 0.1), nrow = 20, ncol = 1)
    colnames(mat) <- "deviance"
    coda::as.mcmc(mat)
  })
  dev_mcmc <- coda::as.mcmc.list(chains)
  result <- occumb:::calc_DIC(dev_mcmc, DIC = TRUE)
  expect_true(!is.null(result))
  expect_equal(names(result), c("pD", "DIC"))
  expect_true(all(is.finite(result)))
})

test_that("calc_DIC returns NULL when deviance contains NA or Inf", {
  na_dev <- make_single_param_mcmc(name = "deviance", fill = NA_real_)
  expect_null(occumb:::calc_DIC(na_dev, DIC = TRUE))
  inf_dev <- make_single_param_mcmc(name = "deviance", fill = Inf)
  expect_null(occumb:::calc_DIC(inf_dev, DIC = TRUE))
})

test_that("order_samples reorders parameters", {
  reordered <- occumb:::order_samples(mcmc2, c("beta", "alpha"))
  expect_equal(
    coda::varnames(reordered),
    c("beta[1]", "beta[2]", "alpha")
  )
})

test_that("order_samples with non-existent params returns empty mcmc.list", {
  result <- occumb:::order_samples(mcmc2, c("nonexistent_only"))
  expect_true(is.null(coda::varnames(result)))
})

test_that("order_samples appends deviance when present but not in params", {
  set.seed(1)
  chains <- lapply(1:2, function(ch) {
    mat <- matrix(rnorm(40), nrow = 20, ncol = 2)
    colnames(mat) <- c("alpha", "deviance")
    coda::as.mcmc(mat)
  })
  mcmc_with_dev <- coda::as.mcmc.list(chains)
  reordered <- occumb:::order_samples(mcmc_with_dev, c("alpha"))
  expect_equal(coda::varnames(reordered), c("alpha", "deviance"))
})

## nimbleSummary() main function -----------------------------------------------

test_that("nimbleSummary returns correct class for mcmc.list input", {
  expect_s3_class(summary2, "nimbleSummary")
  expect_s3_class(summary2, "jagsUI")
})

test_that("nimbleSummary returns correct class for mcmc (single chain) input", {
  set.seed(1)
  mat <- matrix(rnorm(40), nrow = 20, ncol = 2)
  colnames(mat) <- c("a", "b")
  mcmc_obj <- coda::as.mcmc(mat)
  result <- occumb:::nimbleSummary(mcmc_obj)
  expect_s3_class(result, "nimbleSummary")
  expect_equal(result$mcmc.info$n.chains, 1)
})

test_that("nimbleSummary returns correct class for matrix input", {
  set.seed(1)
  mat <- matrix(rnorm(60), nrow = 20, ncol = 3)
  colnames(mat) <- c("a", "b[1]", "b[2]")
  result <- occumb:::nimbleSummary(mat)
  expect_s3_class(result, "nimbleSummary")
  expect_equal(result$mcmc.info$n.chains, 1)
})

test_that("nimbleSummary returns correct class for list input", {
  set.seed(1)
  chain_list <- lapply(1:2, function(ch) {
    mat <- matrix(rnorm(40), nrow = 20, ncol = 2)
    colnames(mat) <- c("x", "y")
    mat
  })
  result <- occumb:::nimbleSummary(chain_list)
  expect_s3_class(result, "nimbleSummary")
  expect_equal(result$mcmc.info$n.chains, 2)
})

test_that("nimbleSummary preserves WAIC from wrapped input", {
  waic_value <- 123.4
  wrapped <- list(samples = mcmc2, WAIC = waic_value)
  result <- occumb:::nimbleSummary(wrapped)
  expect_equal(result$WAIC, waic_value)
})

test_that("nimbleSummary without WAIC has WAIC = NULL", {
  expect_null(summary2$WAIC)
})

test_that("nimbleSummary reorders parameters when specified", {
  result <- occumb:::nimbleSummary(mcmc2, parameters = c("gamma", "alpha"))
  expect_equal(
    rownames(result$summary)[1:3],
    c("gamma[1,1]", "gamma[1,2]", "alpha")
  )
})

test_that("nimbleSummary mcmc.info is correct", {
  expect_equal(summary2$mcmc.info$n.chains, 2)
  expect_equal(summary2$mcmc.info$n.iter, 20)
  expect_equal(summary2$mcmc.info$n.samples, 40)
})

test_that("nimbleSummary errors on unsupported input type", {
  expect_error(occumb:::nimbleSummary("invalid"), "Unsupported input type")
  expect_error(occumb:::nimbleSummary(42), "Unsupported input type")
})

test_that("nimbleSummary output has expected top-level elements", {
  expect_true("sims.list" %in% names(summary2))
  expect_true("summary" %in% names(summary2))
  expect_true("mcmc.info" %in% names(summary2))
  expect_true("samples" %in% names(summary2))
  expect_true(inherits(summary2$samples, "mcmc.list"))
})

## print.nimbleSummary() -------------------------------------------------------

test_that("print.nimbleSummary omits Rhat/n.eff columns for single chain", {
  result <- add_print_fields(summary1)
  output <- capture.output(print(result))
  header_line <- output[grep("mean", output)[1]]
  expect_false(grepl("Rhat", header_line))
})

test_that("print.nimbleSummary includes Rhat/n.eff columns for multiple chains", {
  result <- add_print_fields(summary2)
  output <- capture.output(print(result))
  header_line <- output[grep("mean", output)[1]]
  expect_true(grepl("Rhat", header_line))
  expect_true(grepl("n.eff", header_line))
})

test_that("print.nimbleSummary shows parallel message when parallel = TRUE", {
  result <- add_print_fields(summary2)
  result$parallel <- TRUE
  output <- capture.output(print(result))
  expect_true(any(grepl("ran in parallel", output)))
})

test_that("print.nimbleSummary shows success message when all Rhat < 1.1", {
  result <- add_print_fields(summary2)
  output <- capture.output(print(result))
  expect_true(any(grepl("Successful convergence", output)))
})

test_that("print.nimbleSummary warns about convergence failure", {
  chains <- lapply(1:2, function(ch) {
    mat <- matrix(rnorm(100, mean = ifelse(ch == 1, -100, 100)),
                  nrow = 100, ncol = 1)
    colnames(mat) <- "divergent"
    coda::as.mcmc(mat)
  })
  bad_mcmc <- coda::as.mcmc.list(chains)
  result <- add_print_fields(occumb:::nimbleSummary(bad_mcmc))
  output <- capture.output(print(result))
  expect_true(any(grepl("convergence failure", output)))
})
