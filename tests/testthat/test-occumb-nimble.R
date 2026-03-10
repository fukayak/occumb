### Test cases for write_nimble_model --------------------------------------------
method <- c("pairwise", "full")[1]

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
