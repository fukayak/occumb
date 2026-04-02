### Test cases for write_nimble_model --------------------------------------------
method <- if (identical(Sys.getenv("OCCUMB_TEST_WRITE_NIMBLE_MODEL_FULL"), "true")) "full" else "pairwise"

phi <- theta <- c("i", "ij", "ijk")
psi <- c("i", "ij")
phi_shared <- theta_shared <- psi_shared <- c(FALSE, TRUE)
M_cov_phi <- M_cov_theta <- M_cov_psi <- 1:2
M_cov_phi_shared <- M_cov_theta_shared <- M_cov_psi_shared <- 1:2

factors <- list(
  phi = phi, theta = theta, psi = psi,
  phi_shared = phi_shared, theta_shared = theta_shared, psi_shared = psi_shared,
  M_cov_phi = M_cov_phi, M_cov_theta = M_cov_theta, M_cov_psi = M_cov_psi,
  M_cov_phi_shared = M_cov_phi_shared, M_cov_theta_shared = M_cov_theta_shared,
  M_cov_psi_shared = M_cov_psi_shared
)

if (method == "pairwise") { # 36 cases
  nlevels <- vapply(factors, length, integer(1))
  cases <- suppressMessages(
    DoE.base::oa.design(
      nlevels = nlevels,
      factor.names = factors,
      randomize = FALSE
    )
  )
} else { # 3888 cases
  cases <- do.call(expand.grid, args = factors)
  cases <- subset(cases, phi_shared   | !phi_shared   & M_cov_phi_shared   == 1)
  cases <- subset(cases, theta_shared | !theta_shared & M_cov_theta_shared == 1)
  cases <- subset(cases, psi_shared   | !psi_shared   & M_cov_psi_shared   == 1)
}

### Tests for write_nimble_model() -----------------------------------------------
test_that("NIMBLE code is correct", {
  for (i in seq_len(nrow(cases))) {
    ans <- readLines(system.file("nimble",
                                 "occumb_template1.nimble",
                                 package = "occumb"))
    
    if (cases$phi[i] == "i")
      ans <- c(ans,
               "                r[i, j, k] ~ dgamma(phi[i], 1)")
    if (cases$phi[i] == "ij")
      ans <- c(ans,
               "                r[i, j, k] ~ dgamma(phi[i, j], 1)")
    if (cases$phi[i] == "ijk")
      ans <- c(ans,
               "                r[i, j, k] ~ dgamma(phi[i, j, k], 1)")
    
    ans <- c(ans,
             readLines(system.file("nimble",
                                   "occumb_template2.nimble",
                                   package = "occumb")))
    
    if (cases$theta[i] == "i")
      ans <- c(ans,
               "                u[i, j, k] ~ dbern(z[i, j] * theta[i])")
    if (cases$theta[i] == "ij")
      ans <- c(ans,
               "                u[i, j, k] ~ dbern(z[i, j] * theta[i, j])")
    if (cases$theta[i] == "ijk")
      ans <- c(ans,
               "                u[i, j, k] ~ dbern(z[i, j] * theta[i, j, k])")
    
    ans <- c(ans,
             readLines(system.file("jags",
                                   "occumb_template3.jags",
                                   package = "occumb")))
    
    if (cases$psi[i] == "i")
      ans <- c(ans,
               "            z[i, j] ~ dbern(psi[i])")
    if (cases$psi[i] == "ij")
      ans <- c(ans,
               "            z[i, j] ~ dbern(psi[i, j])")
    
    ans <- c(ans,
             readLines(system.file("jags",
                                   "occumb_template4.jags",
                                   package = "occumb")))
    
    if (cases$phi_shared[i]) {
      if (cases$phi[i] == "i") {
        if (cases$M_cov_phi[i] == 1) {
          term1 <- "alpha[i, 1] * cov_phi[1]"
        } else {
          term1 <- "inprod(alpha[i, 1:M_phi], cov_phi[1:M_phi])"
        }
        if (cases$M_cov_phi_shared[i] == 1) {
          term2 <- "alpha_shared[1] * cov_phi_shared[i, 1]"
        } else {
          term2 <- "inprod(alpha_shared[1:M_phi_shared], cov_phi_shared[i, 1:M_phi_shared])"
        }
        ans <- c(ans, paste0(
          "        log(phi[i]) <- ", term1, " + ", term2))
      } else if (cases$phi[i] == "ij") {
        if (cases$M_cov_phi[i] == 1) {
          term1 <- "alpha[i, 1] * cov_phi[j, 1]"
        } else {
          term1 <- "inprod(alpha[i, 1:M_phi], cov_phi[j, 1:M_phi])"
        }
        if (cases$M_cov_phi_shared[i] == 1) {
          term2 <- "alpha_shared[1] * cov_phi_shared[i, j, 1]"
        } else {
          term2 <- "inprod(alpha_shared[1:M_phi_shared], cov_phi_shared[i, j, 1:M_phi_shared])"
        }
        ans <- c(ans, 
                 "        for (j in 1:J) {", paste0(
                 "            log(phi[i, j]) <- ", term1, " + ", term2),
                 "        }")
      } else if (cases$phi[i] == "ijk") {
        if (cases$M_cov_phi[i] == 1) {
          term1 <- "alpha[i, 1] * cov_phi[j, k, 1]"
        } else {
          term1 <- "inprod(alpha[i, 1:M_phi], cov_phi[j, k, 1:M_phi])"
        }
        if (cases$M_cov_phi_shared[i] == 1) {
          term2 <- "alpha_shared[1] * cov_phi_shared[i, j, k, 1]"
        } else {
          term2 <- "inprod(alpha_shared[1:M_phi_shared], cov_phi_shared[i, j, k, 1:M_phi_shared])"
        }
        ans <- c(ans,
                 "        for (j in 1:J) {",
                 "            for (k in 1:K) {", paste0(
                 "                log(phi[i, j, k]) <- ", term1, " + ", term2),
                 "            }",
                 "        }")
      }
    } else {
      if (cases$phi[i] == "i") {
        if (cases$M_cov_phi[i] == 1) {
          ans <- c(ans,
                   "        log(phi[i]) <- alpha[i, 1] * cov_phi[1]")
        } else {
          ans <- c(ans,
                   "        log(phi[i]) <- inprod(alpha[i, 1:M_phi], cov_phi[1:M_phi])")
        }
      } else if (cases$phi[i] == "ij") {
        if (cases$M_cov_phi[i] == 1) {
          ans <- c(ans,
                   "        for (j in 1:J) {",
                   "            log(phi[i, j]) <- alpha[i, 1] * cov_phi[j, 1]",
                   "        }")
        } else {
          ans <- c(ans,
                   "        for (j in 1:J) {",
                   "            log(phi[i, j]) <- inprod(alpha[i, 1:M_phi], cov_phi[j, 1:M_phi])",
                   "        }")
        }
      } else if (cases$phi[i] == "ijk") {
        if (cases$M_cov_phi[i] == 1) {
          ans <- c(ans,
                   "        for (j in 1:J) {",
                   "            for (k in 1:K) {",
                   "                log(phi[i, j, k]) <- alpha[i, 1] * cov_phi[j, k, 1]",
                   "            }",
                   "        }")
        } else {
          ans <- c(ans,
                   "        for (j in 1:J) {",
                   "            for (k in 1:K) {",
                   "                log(phi[i, j, k]) <- inprod(alpha[i, 1:M_phi], cov_phi[j, k, 1:M_phi])",
                   "            }",
                   "        }")
        }
      }
    }
    
    if (cases$theta_shared[i]) {
      if (cases$theta[i] == "i") {
        if (cases$M_cov_theta[i] == 1) {
          term1 <- "beta[i, 1] * cov_theta[1]"
        } else {
          term1 <- "inprod(beta[i, 1:M_theta], cov_theta[1:M_theta])"
        }
        if (cases$M_cov_theta_shared[i] == 1) {
          term2 <- "beta_shared[1] * cov_theta_shared[i, 1]"
        } else {
          term2 <- "inprod(beta_shared[1:M_theta_shared], cov_theta_shared[i, 1:M_theta_shared])"
        }
        ans <- c(ans, paste0(
          "        logit(theta[i]) <- ", term1, " + ", term2))
      } else if (cases$theta[i] == "ij") {
        if (cases$M_cov_theta[i] == 1) {
          term1 <- "beta[i, 1] * cov_theta[j, 1]"
        } else {
          term1 <- "inprod(beta[i, 1:M_theta], cov_theta[j, 1:M_theta])"
        }
        if (cases$M_cov_theta_shared[i] == 1) {
          term2 <- "beta_shared[1] * cov_theta_shared[i, j, 1]"
        } else {
          term2 <- "inprod(beta_shared[1:M_theta_shared], cov_theta_shared[i, j, 1:M_theta_shared])"
        }
        ans <- c(ans,
                 "        for (j in 1:J) {", paste0(
                 "            logit(theta[i, j]) <- ", term1, " + ", term2),
                 "        }")
      } else if (cases$theta[i] == "ijk") {
        if (cases$M_cov_theta[i] == 1) {
          term1 <- "beta[i, 1] * cov_theta[j, k, 1]"
        } else {
          term1 <- "inprod(beta[i, 1:M_theta], cov_theta[j, k, 1:M_theta])"
        }
        if (cases$M_cov_theta_shared[i] == 1) {
          term2 <- "beta_shared[1] * cov_theta_shared[i, j, k, 1]"
        } else {
          term2 <- "inprod(beta_shared[1:M_theta_shared], cov_theta_shared[i, j, k, 1:M_theta_shared])"
        }
        ans <- c(ans,
                 "        for (j in 1:J) {",
                 "            for (k in 1:K) {", paste0(
                 "                logit(theta[i, j, k]) <- ", term1, " + ", term2),
                 "            }",
                 "        }")
      }
    } else {
      if (cases$theta[i] == "i")
        if (cases$M_cov_theta[i] == 1) {
          ans <- c(ans,
                   "        logit(theta[i]) <- beta[i, 1] * cov_theta[1]")
        } else {
          ans <- c(ans,
                   "        logit(theta[i]) <- inprod(beta[i, 1:M_theta], cov_theta[1:M_theta])")
        }
      if (cases$theta[i] == "ij")
        if (cases$M_cov_theta[i] == 1) {
          ans <- c(ans,
                   "        for (j in 1:J) {",
                   "            logit(theta[i, j]) <- beta[i, 1] * cov_theta[j, 1]",
                   "        }")
        } else {
          ans <- c(ans,
                   "        for (j in 1:J) {",
                   "            logit(theta[i, j]) <- inprod(beta[i, 1:M_theta], cov_theta[j, 1:M_theta])",
                   "        }")
        }
      if (cases$theta[i] == "ijk")
        if (cases$M_cov_theta[i] == 1) {
          ans <- c(ans,
                   "        for (j in 1:J) {",
                   "            for (k in 1:K) {",
                   "                logit(theta[i, j, k]) <- beta[i, 1] * cov_theta[j, k, 1]",
                   "            }",
                   "        }")
        } else {
          ans <- c(ans,
                   "        for (j in 1:J) {",
                   "            for (k in 1:K) {",
                   "                logit(theta[i, j, k]) <- inprod(beta[i, 1:M_theta], cov_theta[j, k, 1:M_theta])",
                   "            }",
                   "        }")
        }
    }
    
    if (cases$psi_shared[i]) {
      if (cases$psi[i] == "i") {
        if (cases$M_cov_psi[i] == 1) {
          term1 <- "gamma[i, 1] * cov_psi[1]"
        } else {
          term1 <- "inprod(gamma[i, 1:M_psi], cov_psi[1:M_psi])"
        }
        if (cases$M_cov_psi_shared[i] == 1) {
          term2 <- "gamma_shared[1] * cov_psi_shared[i, 1]"
        } else {
          term2 <- "inprod(gamma_shared[1:M_psi_shared], cov_psi_shared[i, 1:M_psi_shared])"
        }
        ans <- c(ans,
                 paste0("        logit(psi[i]) <- ", term1, " + ", term2))
        
      } else if (cases$psi[i] == "ij") {
        if (cases$M_cov_psi[i] == 1) {
          term1 <- "gamma[i, 1] * cov_psi[j, 1]"
        } else {
          term1 <- "inprod(gamma[i, 1:M_psi], cov_psi[j, 1:M_psi])"
        }
        if (cases$M_cov_psi_shared[i] == 1) {
          term2 <- "gamma_shared[1] * cov_psi_shared[i, j, 1]"
        } else {
          term2 <- "inprod(gamma_shared[1:M_psi_shared], cov_psi_shared[i, j, 1:M_psi_shared])"
        }
        ans <- c(ans,
                 "        for (j in 1:J) {", paste0(
                 "            logit(psi[i, j]) <- ", term1, " + ", term2),
                 "        }")
      }
    } else {
      if (cases$psi[i] == "i") {
        if (cases$M_cov_psi[i] == 1) {
          ans <- c(ans,
                   "        logit(psi[i]) <- gamma[i, 1] * cov_psi[1]")
        } else {
          ans <- c(ans,
                   "        logit(psi[i]) <- inprod(gamma[i, 1:M_psi], cov_psi[1:M_psi])")
        }
      } else if (cases$psi[i] == "ij") {
        if (cases$M_cov_psi[i] == 1) {
          ans <- c(ans,
                   "        for (j in 1:J) {",
                   "            logit(psi[i, j]) <- gamma[i, 1] * cov_psi[j, 1]",
                   "        }")
        } else {
          ans <- c(ans,
                   "        for (j in 1:J) {",
                   "            logit(psi[i, j]) <- inprod(gamma[i, 1:M_psi], cov_psi[j, 1:M_psi])",
                   "        }")
        }
      }
    }
    
    ans <- c(ans,
             readLines(system.file("nimble",
                                   "occumb_template5.nimble",
                                   package = "occumb")))
    
    if (cases$phi_shared[i])
      ans <- c(ans,
               "    for (m in 1:M_phi_shared) {",
               "        alpha_shared[m] ~ dnorm(0, prior_prec)",
               "    }")
    if (cases$theta_shared[i])
      ans <- c(ans,
               "    for (m in 1:M_theta_shared) {",
               "        beta_shared[m] ~ dnorm(0, prior_prec)",
               "    }")
    if (cases$psi_shared[i])
      ans <- c(ans,
               "    for (m in 1:M_psi_shared) {",
               "        gamma_shared[m] ~ dnorm(0, prior_prec)",
               "    }")
    
    ans <- c(ans, "}", "")
    
    res <- write_nimble_model(phi                = cases$phi[i],
                              theta              = cases$theta[i],
                              psi                = cases$psi[i],
                              phi_shared         = cases$phi_shared[i],
                              theta_shared       = cases$theta_shared[i],
                              psi_shared         = cases$psi_shared[i],
                              M_cov_phi          = cases$M_cov_phi[i], 
                              M_cov_phi_shared   = cases$M_cov_phi_shared[i],
                              M_cov_theta        = cases$M_cov_theta[i],
                              M_cov_theta_shared = cases$M_cov_theta_shared[i],
                              M_cov_psi          = cases$M_cov_psi[i],
                              M_cov_psi_shared   = cases$M_cov_psi_shared[i])
    expect_equal(res, ans)
  }
})

### Tests for NIMBLE setup helper functions ------------------------------------

make_test_const <- function(I = 2, J = 2, K = 2) {
  list(I = I, J = J, K = K,
       N = matrix(10L, nrow = J, ncol = K))
}

make_test_data <- function(m_phi = 1L, m_theta = 2L, m_psi = 3L,
                           shared = FALSE) {
  I <- 2; J <- 2; K <- 2
  M <- max(c(m_phi, m_theta, m_psi))
  data <- list(
    y = array(1L, dim = c(I, J, K)),
    M = M,
    m_phi = m_phi, m_theta = m_theta, m_psi = m_psi,
    prior_prec = 0.01, prior_ulim = 20,
    cov_phi = matrix(1, J, 1),
    cov_theta = matrix(1, J, 1),
    cov_psi = matrix(1, J, 1)
  )
  if (shared) {
    data$M_phi_shared <- 2L
    data$M_theta_shared <- 2L
    data$M_psi_shared <- 2L
    data$cov_phi_shared <- matrix(1, I, 2)
    data$cov_theta_shared <- matrix(1, J, 2)
    data$cov_psi_shared <- matrix(1, I, 2)
  }
  data
}

default_const <- make_test_const()
default_data <- make_test_data()
default_data_shared <- make_test_data(shared = TRUE)
default_const_result <- occumb:::set_const_nimble(default_const, default_data)
default_data_result <- occumb:::set_data_nimble(default_data)
mock_inits <- function() list(z = 1, rho = matrix(0, 1, 2))

## pad_dummy_value / make_rho_index --------------------------------------------

test_that("pad_dummy_value pads length-1 vector with -999", {
  expect_equal(occumb:::pad_dummy_value(5L), c(5L, -999))
})

test_that("pad_dummy_value leaves length-2+ vector unchanged", {
  expect_equal(occumb:::pad_dummy_value(c(1L, 2L)), c(1L, 2L))
  expect_equal(occumb:::pad_dummy_value(1:5), 1:5)
})

test_that("make_rho_index produces correct upper-triangular index matrix", {
  idx <- occumb:::make_rho_index(3)
  expected <- matrix(c(0L, 0L, 0L,
                       1L, 0L, 0L,
                       2L, 3L, 0L), nrow = 3, ncol = 3)
  expect_equal(idx, expected)
})

## set_const_nimble() ----------------------------------------------------------

test_that("set_const_nimble extracts I, J, K, N, M, prior_prec, prior_ulim", {
  const <- make_test_const(I = 3, J = 4, K = 5)
  result <- occumb:::set_const_nimble(const, default_data)
  expect_equal(result$I, 3)
  expect_equal(result$J, 4)
  expect_equal(result$K, 5)
  expect_equal(result$N, matrix(10L, nrow = 4, ncol = 5))
  expect_equal(result$M, default_data$M)
  expect_equal(result$prior_prec, default_data$prior_prec)
  expect_equal(result$prior_ulim, default_data$prior_ulim)
})

test_that("set_const_nimble computes M_phi/M_theta/M_psi from m_* lengths", {
  data <- make_test_data(m_phi = 1L, m_theta = 2:3, m_psi = 4:6)
  result <- occumb:::set_const_nimble(default_const, data)
  expect_equal(result$M_phi, 1)
  expect_equal(result$M_theta, 2)
  expect_equal(result$M_psi, 3)
})

test_that("set_const_nimble pads length-1 m_* with dummy value", {
  expect_equal(default_const_result$m_phi, c(1L, -999))
  expect_equal(default_const_result$m_theta, c(2L, -999))
  expect_equal(default_const_result$m_psi, c(3L, -999))
})

test_that("set_const_nimble does not pad length-2+ m_*", {
  data <- make_test_data(m_phi = 1:2, m_theta = 3:4, m_psi = 5:6)
  result <- occumb:::set_const_nimble(default_const, data)
  expect_equal(result$m_phi, 1:2)
  expect_equal(result$m_theta, 3:4)
  expect_equal(result$m_psi, 5:6)
})

test_that("set_const_nimble includes correct rho_index", {
  expect_equal(default_const_result$rho_index,
               occumb:::make_rho_index(default_data$M))
})

test_that("set_const_nimble excludes M_*_shared when not shared", {
  expect_false("M_phi_shared" %in% names(default_const_result))
  expect_false("M_theta_shared" %in% names(default_const_result))
  expect_false("M_psi_shared" %in% names(default_const_result))
})

test_that("set_const_nimble includes M_*_shared when shared", {
  result <- occumb:::set_const_nimble(default_const, default_data_shared)
  expect_equal(result$M_phi_shared, 2L)
  expect_equal(result$M_theta_shared, 2L)
  expect_equal(result$M_psi_shared, 2L)
})

## set_data_nimble() -----------------------------------------------------------

test_that("set_data_nimble extracts y and cov_* with correct values", {
  expect_equal(default_data_result$y, default_data$y)
  expect_equal(default_data_result$cov_phi, default_data$cov_phi)
  expect_equal(default_data_result$cov_theta, default_data$cov_theta)
  expect_equal(default_data_result$cov_psi, default_data$cov_psi)
})

test_that("set_data_nimble excludes dimension and prior elements", {
  excluded <- c("M", "m_phi", "m_theta", "m_psi", "prior_prec", "prior_ulim")
  for (nm in excluded) {
    expect_false(nm %in% names(default_data_result),
                 info = paste(nm, "should be excluded"))
  }
})

test_that("set_data_nimble excludes cov_*_shared when not shared", {
  expect_false("cov_phi_shared" %in% names(default_data_result))
  expect_false("cov_theta_shared" %in% names(default_data_result))
  expect_false("cov_psi_shared" %in% names(default_data_result))
})

test_that("set_data_nimble includes cov_*_shared when shared", {
  result <- occumb:::set_data_nimble(default_data_shared)
  expect_true("cov_phi_shared" %in% names(result))
  expect_true("cov_theta_shared" %in% names(result))
  expect_true("cov_psi_shared" %in% names(result))
})

## set_inits_nimble() ----------------------------------------------------------

test_that("set_inits_nimble with seed=TRUE produces reproducible results", {
  skip_if_not_installed("nimble")
  result1 <- occumb:::set_inits_nimble(mock_inits, seed = TRUE,
                                       n.chains = 2, n_rho = 3)
  result2 <- occumb:::set_inits_nimble(mock_inits, seed = TRUE,
                                       n.chains = 2, n_rho = 3)
  expect_equal(result1[[1]]$z, result2[[1]]$z)
  expect_equal(result1[[2]]$z, result2[[2]]$z)
  expect_equal(result1[[1]]$rho, result2[[1]]$rho)
  expect_equal(result1[[1]]$.RNG.seed, result2[[1]]$.RNG.seed)
  expect_equal(result1[[2]]$.RNG.seed, result2[[2]]$.RNG.seed)
})

test_that("set_inits_nimble with seed=NULL returns n.chains elements", {
  skip_if_not_installed("nimble")
  result <- occumb:::set_inits_nimble(mock_inits, seed = NULL,
                                      n.chains = 3, n_rho = 1)
  expect_length(result, 3)
})

test_that("set_inits_nimble with seed=FALSE produces distinct per-chain seeds", {
  skip_if_not_installed("nimble")
  result <- occumb:::set_inits_nimble(mock_inits, seed = FALSE,
                                      n.chains = 3, n_rho = 1)
  expect_length(result, 3)
  seeds <- vapply(result, function(x) x$.RNG.seed, numeric(1))
  expect_equal(length(unique(seeds)), 3)
})

test_that("set_inits_nimble with numeric vector uses seeds directly", {
  skip_if_not_installed("nimble")
  result <- occumb:::set_inits_nimble(mock_inits, seed = c(10, 20),
                                      n.chains = 2, n_rho = 1)
  expect_equal(result[[1]]$.RNG.seed, 10)
  expect_equal(result[[2]]$.RNG.seed, 20)
})

test_that("set_inits_nimble with single numeric expands to per-chain seeds", {
  skip_if_not_installed("nimble")
  expect_message(
    result <- occumb:::set_inits_nimble(mock_inits, seed = 100,
                                        n.chains = 3, n_rho = 1),
    "Expanding"
  )
  expect_equal(result[[1]]$.RNG.seed, 100)
  expect_equal(result[[2]]$.RNG.seed, 101)
  expect_equal(result[[3]]$.RNG.seed, 102)
})

test_that("set_inits_nimble errors on invalid seed", {
  skip_if_not_installed("nimble")
  expect_error(
    occumb:::set_inits_nimble(mock_inits, seed = "abc",
                              n.chains = 2, n_rho = 1),
    "Invalid 'seed'"
  )
})

test_that("set_inits_nimble errors when seed vector length != n.chains", {
  skip_if_not_installed("nimble")
  expect_error(
    occumb:::set_inits_nimble(mock_inits, seed = c(1, 2, 3),
                              n.chains = 2, n_rho = 1),
    "Invalid 'seed'"
  )
})

test_that("set_inits_nimble replaces rho with double(n_rho)", {
  skip_if_not_installed("nimble")
  mock_inits_large_rho <- function() list(z = 1, rho = matrix(99, 2, 3))
  result <- occumb:::set_inits_nimble(mock_inits_large_rho, seed = TRUE,
                                      n.chains = 1, n_rho = 5)
  expect_equal(result[[1]]$rho, double(5))
})

test_that("set_inits_nimble includes .RNG.seed", {
  skip_if_not_installed("nimble")
  result <- occumb:::set_inits_nimble(mock_inits, seed = TRUE,
                                      n.chains = 2, n_rho = 1)
  for (i in seq_along(result)) {
    expect_true(".RNG.seed" %in% names(result[[i]]))
    expect_true(is.numeric(result[[i]]$.RNG.seed))
  }
})

test_that("set_inits_nimble warns on duplicate seeds", {
  skip_if_not_installed("nimble")
  old_verbose <- nimble::getNimbleOption("verbose")
  nimble::nimbleOptions(verbose = TRUE)
  on.exit(nimble::nimbleOptions(verbose = old_verbose))
  expect_message(
    occumb:::set_inits_nimble(mock_inits, seed = c(42, 42),
                              n.chains = 2, n_rho = 1),
    "duplicates"
  )
})
