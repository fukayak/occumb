### Tests for set_const() ------------------------------------------------------
test_that("Replacement of missing y to NA works", {
    y         <- array(1:8, dim = rep(2, 3))
    y[, 1, 1] <- 0
    result    <- set_const(occumbData(y = y))
    expect_identical(result$y[, 1, 1], as.numeric(rep(NA, 2)))
    expect_identical(result$y[, 2, 1], c(3, 4))
    expect_identical(result$y[, 1, 2], c(5, 6))
    expect_identical(result$y[, 2, 2], c(7, 8))
})
test_that("Replacement of sequence depth 0 to 1 works", {
    y         <- array(1:8, dim = rep(2, 3))
    y[, 1, 1] <- 0
    result    <- set_const(occumbData(y = y))
    expect_equal(result$N[1, 1], 1)
    expect_equal(result$N[2, 1], sum(3:4))
    expect_equal(result$N[1, 2], sum(5:6))
    expect_equal(result$N[2, 2], sum(7:8))
})

### Tests for set_modargs() ----------------------------------------------------
I <- 2; J <- 3; K <- 4
cov1 <- rnorm(I)
cov2 <- factor(1:I)
cov3 <- rnorm(J)
cov4 <- factor(1:J)
cov5 <- matrix(rnorm(J * K), nrow = J)
cov6 <- matrix(rep(factor(1:K), each = J), nrow = J)
data <- occumbData(y = array(0, dim = c(I, J, K)),
                   spec_cov = list(cov1 = cov1, cov2 = cov2),
                   site_cov = list(cov3 = cov3, cov4 = cov4),
                   repl_cov = list(cov5 = cov5, cov6 = cov6))

test_that("Temp: psi correct", {
    ## Output
    # null model
    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, NULL, data)
    expect_equal(result$psi, "i")
    expect_equal(result$M, 3)
    expect_equal(result$cov_psi, 1)
    expect_equal(result$m_psi, 3)

    # spec_cov (not allowed)
    expect_error(set_modargs(~ 1, ~ 1, ~ cov1, NULL, NULL, NULL, data),
                 sprintf("Unexpected terms in formula_psi: %s
Note that only site covariates are allowed for formula_psi.",
                         "cov1"))

    # site_cov (continuous)
    result <- set_modargs(~ 1, ~ 1, ~ cov3, NULL, NULL, NULL, data)
    expect_equal(result$psi, "ij")
    expect_equal(result$M, 4)
    ans_cov <- array(dim = c(I, J, 2))
    for (i in 1:I) {
        for (j in 1:J) {
        ans_cov[i, j, 1] <- 1
        ans_cov[i, j, 2] <- cov3[j]
        }
    }
    dimnames(ans_cov)[[3]] <- c("(Intercept)", "cov3")
    expect_equal(result$cov_psi, ans_cov)
    expect_equal(result$m_psi, 3:4)

    # site_cov (factor)
    result <- set_modargs(~ 1, ~ 1, ~ cov4, NULL, NULL, NULL, data)
    expect_equal(result$psi, "ij")
    expect_equal(result$M, 5)
    ans_cov <- array(dim = c(I, J, 3))
    for (i in 1:I) {
        for (j in 1:J) {
        ans_cov[i, j, 1] <- 1
        ans_cov[i, j, 2] <- as.numeric(cov4[j] == 2)
        ans_cov[i, j, 3] <- as.numeric(cov4[j] == 3)
        }
    }
    dimnames(ans_cov)[[3]] <- c("(Intercept)", "cov42", "cov43")
    expect_equal(result$cov_psi, ans_cov)
    expect_equal(result$m_psi, 3:5)

    # site_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ cov3 * cov4, NULL, NULL, NULL, data)
    expect_equal(result$psi, "ij")
    expect_equal(result$M, 8)
    ans_cov <- array(dim = c(I, J, 6))
    for (i in 1:I) {
        for (j in 1:J) {
        ans_cov[i, j, 1] <- 1
        ans_cov[i, j, 2] <- cov3[j]
        ans_cov[i, j, 3] <- as.numeric(cov4[j] == 2)
        ans_cov[i, j, 4] <- as.numeric(cov4[j] == 3)
        ans_cov[i, j, 5] <- cov3[j] * as.numeric(cov4[j] == 2)
        ans_cov[i, j, 6] <- cov3[j] * as.numeric(cov4[j] == 3)
        }
    }
    dimnames(ans_cov)[[3]] <- c("(Intercept)", "cov3", "cov42", "cov43",
                                "cov3:cov42", "cov3:cov43")
    expect_equal(result$cov_psi, ans_cov)
    expect_equal(result$m_psi, 3:8)

    # repl_cov (not allowed)
    expect_error(set_modargs(~ 1, ~ 1, ~ cov5, NULL, NULL, NULL, data),
                 sprintf("Unexpected terms in formula_psi: %s
Note that only site covariates are allowed for formula_psi.",
                         "cov5"))

    ## Errors and Warnings
    expect_error(set_modargs(~ 1, ~ 1, ~ 0, NULL, NULL, NULL, data),
                 sprintf("No intercept in formula_%s: remove 0 or -1 from the formula", "psi"))
    expect_error(set_modargs(~ 1, ~ 1, ~ -1, NULL, NULL, NULL, data),
                 sprintf("No intercept in formula_%s: remove 0 or -1 from the formula", "psi"))
    expect_error(set_modargs(~ 1, ~ 1, ~ xxx, NULL, NULL, NULL, data),
                 sprintf("Unexpected terms in formula_psi: %s
Note that only site covariates are allowed for formula_psi.",
                         "xxx"))
})
test_that("Temp: psi_shared correct", {
    ## Output
    # spec_cov (continuous)
    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov1, data)
    expect_equal(result$psi, "i")
    expect_true(result$psi_shared)
    expect_equal(result$M_psi_shared, 1)
    ans_cov <- array(dim = c(I, 1))
    for (i in 1:I) {
        ans_cov[i, 1] <- cov1[i]
    }
    colnames(ans_cov) <- "cov1"
    expect_equal(result$cov_psi_shared, ans_cov)

    # spec_cov (factor)
    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov2, data)
    expect_equal(result$psi, "i")
    expect_true(result$psi_shared)
    expect_equal(result$M_psi_shared, 1)
    ans_cov <- array(dim = c(I, 1))
    for (i in 1:I) {
        ans_cov[i, 1] <- as.numeric((cov2[i] == 2))
    }
    colnames(ans_cov) <- c("cov22")
    expect_equal(result$cov_psi_shared, ans_cov)

    # spec_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov1 * cov2, data)
    expect_equal(result$psi, "i")
    expect_true(result$psi_shared)
    expect_equal(result$M_psi_shared, 3)
    ans_cov <- array(dim = c(I, 3))
    for (i in 1:I) {
        ans_cov[i, 1] <- cov1[i]
        ans_cov[i, 2] <- as.numeric((cov2[i] == 2))
        ans_cov[i, 3] <- cov1[i] * as.numeric((cov2[i] == 2))
    }
    colnames(ans_cov) <- c("cov1", "cov22", "cov1:cov22")
    expect_equal(result$cov_psi_shared, ans_cov)

    # site_cov (continuous)
    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov3, data)
    expect_equal(result$psi, "ij")
    expect_true(result$psi_shared)
    expect_equal(result$M_psi_shared, 1)
    ans_cov <- array(dim = c(I, J, 1))
    for (i in 1:I) {
        for (j in 1:J) {
        ans_cov[i, j, 1] <- cov3[j]
        }
    }
    dimnames(ans_cov)[[3]] <- "cov3"
    expect_equal(result$cov_psi_shared, ans_cov)

    # site_cov (factor)
    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov4, data)
    expect_true(result$psi_shared)
    expect_equal(result$M_psi_shared, 2)
    ans_cov <- array(dim = c(I, J, 2))
    for (i in 1:I) {
        for (j in 1:J) {
        ans_cov[i, j, 1] <- as.numeric(cov4[j] == 2)
        ans_cov[i, j, 2] <- as.numeric(cov4[j] == 3)
        }
    }
    dimnames(ans_cov)[[3]] <- c("cov42", "cov43")
    expect_equal(result$cov_psi_shared, ans_cov)

    # site_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov3 * cov4, data)
    expect_true(result$psi_shared)
    expect_equal(result$M_psi_shared, 5)
    ans_cov <- array(dim = c(I, J, 5))
    for (i in 1:I) {
        for (j in 1:J) {
        ans_cov[i, j, 1] <- cov3[j]
        ans_cov[i, j, 2] <- as.numeric(cov4[j] == 2)
        ans_cov[i, j, 3] <- as.numeric(cov4[j] == 3)
        ans_cov[i, j, 4] <- cov3[j] * as.numeric(cov4[j] == 2)
        ans_cov[i, j, 5] <- cov3[j] * as.numeric(cov4[j] == 3)
        }
    }
    dimnames(ans_cov)[[3]] <- c("cov3", "cov42", "cov43",
                                "cov3:cov42", "cov3:cov43")
    expect_equal(result$cov_psi_shared, ans_cov)

    # spec_cov * site_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov1 * cov4, data)
    expect_equal(result$psi, "ij")
    expect_true(result$psi_shared)
    expect_equal(result$M_psi_shared, 5)
    ans_cov <- array(dim = c(I, J, 5))
    for (i in 1:I) {
        for (j in 1:J) {
            ans_cov[i, j, 1] <- cov1[i]
            ans_cov[i, j, 2] <- as.numeric(cov4[j] == 2)
            ans_cov[i, j, 3] <- as.numeric(cov4[j] == 3)
            ans_cov[i, j, 4] <- cov1[i] * as.numeric(cov4[j] == 2)
            ans_cov[i, j, 5] <- cov1[i] * as.numeric(cov4[j] == 3)
        }
    }
    dimnames(ans_cov)[[3]] <- c("cov1", "cov42", "cov43",
                                "cov1:cov42", "cov1:cov43")
    expect_equal(result$cov_psi_shared, ans_cov)

    # repl_cov (not allowed)
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov5, data),
                             sprintf("Unexpected terms in formula_psi_shared: %s
Note that only site covariates, species covariates, or their interactions are allowed for formula_psi_shared.", "cov5"))

    ## Errors and Warnings
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ 0, data),
                 sprintf("No intercept in formula_%s: remove 0 or -1 from the formula", "psi_shared"))
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ -1, data),
                 sprintf("No intercept in formula_%s: remove 0 or -1 from the formula", "psi_shared"))
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ xxx, data),
                             sprintf("Unexpected terms in formula_psi_shared: %s
Note that only site covariates, species covariates, or their interactions are allowed for formula_psi_shared.", "xxx"))
})

test_that("Temp: theta correct", {
    ## Output
    # null model
    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, NULL, data)
    expect_equal(result$theta, "i")
    expect_equal(result$M, 3)
    expect_equal(result$cov_theta, 1)
    expect_equal(result$m_theta, 2)

    # spec_cov (not allowed)
    expect_error(set_modargs(~ 1, ~ cov1, ~ 1, NULL, NULL, NULL, data),
                 sprintf("Unexpected terms in formula_theta: %s
Note that species covariates are not allowed for formula_theta.",
                         "cov1"))

    # site_cov (continuous)
    result <- set_modargs(~ 1, ~ cov3, ~ 1, NULL, NULL, NULL, data)
    expect_equal(result$theta, "ij")
    expect_equal(result$M, 4)
    ans_cov <- array(dim = c(I, J, 2))
    for (i in 1:I) {
        for (j in 1:J) {
        ans_cov[i, j, 1] <- 1
        ans_cov[i, j, 2] <- cov3[j]
        }
    }
    dimnames(ans_cov)[[3]] <- c("(Intercept)", "cov3")
    expect_equal(result$cov_theta, ans_cov)
    expect_equal(result$m_theta, 2:3)

    # site_cov (factor)
    result <- set_modargs(~ 1, ~ cov4, ~ 1, NULL, NULL, NULL, data)
    expect_equal(result$theta, "ij")
    expect_equal(result$M, 5)
    ans_cov <- array(dim = c(I, J, 3))
    for (i in 1:I) {
        for (j in 1:J) {
        ans_cov[i, j, 1] <- 1
        ans_cov[i, j, 2] <- as.numeric((cov4[j] == 2))
        ans_cov[i, j, 3] <- as.numeric((cov4[j] == 3))
        }
    }
    dimnames(ans_cov)[[3]] <- c("(Intercept)", "cov42", "cov43")
    expect_equal(result$cov_theta, ans_cov)
    expect_equal(result$m_theta, 2:4)

    # site_cov (interaction)
    result <- set_modargs(~ 1, ~ cov3 * cov4, ~ 1, NULL, NULL, NULL, data)
    expect_equal(result$theta, "ij")
    expect_equal(result$M, 8)
    ans_cov <- array(dim = c(I, J, 6))
    for (i in 1:I) {
        for (j in 1:J) {
        ans_cov[i, j, 1] <- 1
        ans_cov[i, j, 2] <- cov3[j]
        ans_cov[i, j, 3] <- as.numeric(cov4[j] == 2)
        ans_cov[i, j, 4] <- as.numeric(cov4[j] == 3)
        ans_cov[i, j, 5] <- cov3[j] * as.numeric(cov4[j] == 2)
        ans_cov[i, j, 6] <- cov3[j] * as.numeric(cov4[j] == 3)
        }
    }
    dimnames(ans_cov)[[3]] <- c("(Intercept)", "cov3", "cov42", "cov43",
                                "cov3:cov42", "cov3:cov43")
    expect_equal(result$cov_theta, ans_cov)
    expect_equal(result$m_theta, 2:7)

    # repl_cov (continuous)
    result <- set_modargs(~ 1, ~ cov5, ~ 1, NULL, NULL, NULL, data)
    expect_equal(result$theta, "ijk")
    expect_equal(result$M, 4)
    ans_cov <- array(dim = c(I, J, K, 2))
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                ans_cov[i, j, k, 1] <- 1
                ans_cov[i, j, k, 2] <- cov5[j, k]
            }
        }
    }
    dimnames(ans_cov)[[4]] <- c("(Intercept)", "cov5")
    expect_equal(result$cov_theta, ans_cov)
    expect_equal(result$m_theta, 2:3)

    # repl_cov (factor)
    result <- set_modargs(~ 1, ~ cov6, ~ 1, NULL, NULL, NULL, data)
    expect_equal(result$theta, "ijk")
    expect_equal(result$M, 6)
    ans_cov <- array(dim = c(I, J, K, 4))
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                ans_cov[i, j, k, 1] <- 1
                ans_cov[i, j, k, 2] <- as.numeric((cov6[j, k] == 2))
                ans_cov[i, j, k, 3] <- as.numeric((cov6[j, k] == 3))
                ans_cov[i, j, k, 4] <- as.numeric((cov6[j, k] == 4))
            }
        }
    }
    dimnames(ans_cov)[[4]] <- c("(Intercept)", "cov62", "cov63", "cov64")
    expect_equal(result$cov_theta, ans_cov)
    expect_equal(result$m_theta, 2:5)

    # repl_cov (interaction)
    result <- set_modargs(~ 1, ~ cov5 * cov6, ~ 1, NULL, NULL, NULL, data)
    expect_equal(result$theta, "ijk")
    expect_equal(result$M, 10)
    ans_cov <- array(dim = c(I, J, K, 8))
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                ans_cov[i, j, k, 1] <- 1
                ans_cov[i, j, k, 2] <- cov5[j, k]
                ans_cov[i, j, k, 3] <- as.numeric(cov6[j, k] == 2)
                ans_cov[i, j, k, 4] <- as.numeric(cov6[j, k] == 3)
                ans_cov[i, j, k, 5] <- as.numeric(cov6[j, k] == 4)
                ans_cov[i, j, k, 6] <- cov5[j, k] * as.numeric(cov6[j, k] == 2)
                ans_cov[i, j, k, 7] <- cov5[j, k] * as.numeric(cov6[j, k] == 3)
                ans_cov[i, j, k, 8] <- cov5[j, k] * as.numeric(cov6[j, k] == 4)
            }
        }
    }
    dimnames(ans_cov)[[4]] <- c("(Intercept)", "cov5", "cov62", "cov63", "cov64",
                                "cov5:cov62", "cov5:cov63", "cov5:cov64")
    expect_equal(result$cov_theta, ans_cov)
    expect_equal(result$m_theta, 2:9)

    # site_cov * repl_cov (interaction)
    result <- set_modargs(~ 1, ~ cov3 * cov6, ~ 1, NULL, NULL, NULL, data)
    expect_equal(result$theta, "ijk")
    expect_equal(result$M, 10)
    ans_cov <- array(dim = c(I, J, K, 8))
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                ans_cov[i, j, k, 1] <- 1
                ans_cov[i, j, k, 2] <- cov3[j]
                ans_cov[i, j, k, 3] <- as.numeric(cov6[j, k] == 2)
                ans_cov[i, j, k, 4] <- as.numeric(cov6[j, k] == 3)
                ans_cov[i, j, k, 5] <- as.numeric(cov6[j, k] == 4)
                ans_cov[i, j, k, 6] <- cov3[j] * as.numeric(cov6[j, k] == 2)
                ans_cov[i, j, k, 7] <- cov3[j] * as.numeric(cov6[j, k] == 3)
                ans_cov[i, j, k, 8] <- cov3[j] * as.numeric(cov6[j, k] == 4)
            }
        }
    }
    dimnames(ans_cov)[[4]] <- c("(Intercept)", "cov3", "cov62", "cov63", "cov64",
                                "cov3:cov62", "cov3:cov63", "cov3:cov64")
    expect_equal(result$cov_theta, ans_cov)
    expect_equal(result$m_theta, 2:9)

    ## Errors and Warnings
    expect_error(set_modargs(~ 1, ~ 0, ~ 1, NULL, NULL, NULL, data),
                 sprintf("No intercept in formula_%s: remove 0 or -1 from the formula", "theta"))
    expect_error(set_modargs(~ 1, ~ -1, ~ 1, NULL, NULL, NULL, data),
                 sprintf("No intercept in formula_%s: remove 0 or -1 from the formula", "theta"))
    expect_error(set_modargs(~ 1, ~ xxx, ~ 1, NULL, NULL, NULL, data),
                 sprintf("Unexpected terms in formula_theta: %s
Note that species covariates are not allowed for formula_theta.",
                         "xxx"))
})
test_that("Temp: theta_shared correct", {
#    ## Output
#    # spec_cov (continuous)
#    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov1, data)
#    expect_equal(result$psi, "i")
#    expect_true(result$psi_shared)
#    expect_equal(result$M_psi_shared, 1)
#    ans_cov <- array(dim = c(I, 1))
#    for (i in 1:I) {
#        ans_cov[i, 1] <- cov1[i]
#    }
#    colnames(ans_cov) <- "cov1"
#    expect_equal(result$cov_psi_shared, ans_cov)
#
#    # spec_cov (factor)
#    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov2, data)
#    expect_equal(result$psi, "i")
#    expect_true(result$psi_shared)
#    expect_equal(result$M_psi_shared, 1)
#    ans_cov <- array(dim = c(I, 1))
#    for (i in 1:I) {
#        ans_cov[i, 1] <- as.numeric((cov2[i] == 2))
#    }
#    colnames(ans_cov) <- c("cov22")
#    expect_equal(result$cov_psi_shared, ans_cov)
#
#    # spec_cov (interaction)
#    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov1 * cov2, data)
#    expect_equal(result$psi, "i")
#    expect_true(result$psi_shared)
#    expect_equal(result$M_psi_shared, 3)
#    ans_cov <- array(dim = c(I, 3))
#    for (i in 1:I) {
#        ans_cov[i, 1] <- cov1[i]
#        ans_cov[i, 2] <- as.numeric((cov2[i] == 2))
#        ans_cov[i, 3] <- cov1[i] * as.numeric((cov2[i] == 2))
#    }
#    colnames(ans_cov) <- c("cov1", "cov22", "cov1:cov22")
#    expect_equal(result$cov_psi_shared, ans_cov)
#
#    # site_cov (continuous)
#    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov3, data)
#    expect_equal(result$psi, "ij")
#    expect_true(result$psi_shared)
#    expect_equal(result$M_psi_shared, 1)
#    ans_cov <- array(dim = c(I, J, 1))
#    for (i in 1:I) {
#        for (j in 1:J) {
#        ans_cov[i, j, 1] <- cov3[j]
#        }
#    }
#    dimnames(ans_cov)[[3]] <- "cov3"
#    expect_equal(result$cov_psi_shared, ans_cov)
#
#    # site_cov (factor)
#    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov4, data)
#    expect_true(result$psi_shared)
#    expect_equal(result$M_psi_shared, 2)
#    ans_cov <- array(dim = c(I, J, 2))
#    for (i in 1:I) {
#        for (j in 1:J) {
#        ans_cov[i, j, 1] <- as.numeric(cov4[j] == 2)
#        ans_cov[i, j, 2] <- as.numeric(cov4[j] == 3)
#        }
#    }
#    dimnames(ans_cov)[[3]] <- c("cov42", "cov43")
#    expect_equal(result$cov_psi_shared, ans_cov)
#
#    # site_cov (interaction)
#    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov3 * cov4, data)
#    expect_true(result$psi_shared)
#    expect_equal(result$M_psi_shared, 5)
#    ans_cov <- array(dim = c(I, J, 5))
#    for (i in 1:I) {
#        for (j in 1:J) {
#        ans_cov[i, j, 1] <- cov3[j]
#        ans_cov[i, j, 2] <- as.numeric(cov4[j] == 2)
#        ans_cov[i, j, 3] <- as.numeric(cov4[j] == 3)
#        ans_cov[i, j, 4] <- cov3[j] * as.numeric(cov4[j] == 2)
#        ans_cov[i, j, 5] <- cov3[j] * as.numeric(cov4[j] == 3)
#        }
#    }
#    dimnames(ans_cov)[[3]] <- c("cov3", "cov42", "cov43",
#                                "cov3:cov42", "cov3:cov43")
#    expect_equal(result$cov_psi_shared, ans_cov)
#
#    # spec_cov * site_cov (interaction)
#    result <- set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov1 * cov4, data)
#    expect_equal(result$psi, "ij")
#    expect_true(result$psi_shared)
#    expect_equal(result$M_psi_shared, 5)
#    ans_cov <- array(dim = c(I, J, 5))
#    for (i in 1:I) {
#        for (j in 1:J) {
#            ans_cov[i, j, 1] <- cov1[i]
#            ans_cov[i, j, 2] <- as.numeric(cov4[j] == 2)
#            ans_cov[i, j, 3] <- as.numeric(cov4[j] == 3)
#            ans_cov[i, j, 4] <- cov1[i] * as.numeric(cov4[j] == 2)
#            ans_cov[i, j, 5] <- cov1[i] * as.numeric(cov4[j] == 3)
#        }
#    }
#    dimnames(ans_cov)[[3]] <- c("cov1", "cov42", "cov43",
#                                "cov1:cov42", "cov1:cov43")
#    expect_equal(result$cov_psi_shared, ans_cov)
#
#    # repl_cov (not allowed)
#    expect_error(set_modargs(~ 1, ~ 1, ~ 1, NULL, NULL, ~ cov5, data),
#                             sprintf("Unexpected terms in formula_psi_shared: %s
#Note that only site covariates, species covariates, or their interactions are allowed for formula_psi_shared.", "cov5"))

    ## Errors and Warnings
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, NULL, ~ 0, NULL, data),
                 "No intercept in formula_theta_shared: remove 0 or -1 from the formula")
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, NULL, ~ -1, NULL, data),
                 "No intercept in formula_theta_shared: remove 0 or -1 from the formula")
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, NULL, ~ xxx, NULL, data),
                             sprintf("Unexpected terms in formula_theta_shared: %s
Make sure they are found in either spec_cov, site_cov, or repl_cov.", "xxx"))
})


#test_that("Setup for a null model works", {
#    result <- set_modargs(~ 1, ~ 1, ~ 1,
#                          occumbData(y = array(0, dim = rep(2, 3))))
#    expect_equal(result$M, 3)
#    expect_equal(result$code, "111000")
#    expect_equal(result$cov_phi, 1)
#    expect_equal(result$cov_theta, 1)
#    expect_equal(result$cov_psi, 1)
#    expect_equal(result$m_phi, 1)
#    expect_equal(result$m_theta, 2)
#    expect_equal(result$m_psi, 3)
#})
#test_that("Setup for varying-psi model works", {
#    cov1 <- rnorm(2); cov2 <- factor(1:2)
#    data <- occumbData(y = array(0, dim = rep(2, 3)),
#                       site_cov = list(cov1 = cov1,
#                                       cov2 = cov2))
#
#    # Continuous effect
#    result <- set_modargs(~ 1, ~ 1, ~ cov1, data)
#    mm     <- model.matrix(as.formula(~ cov1), data@site_cov)
#    expect_equal(result$M, 4)
#    expect_equal(result$code, "112000")
#    expect_equal(result$cov_phi, 1)
#    expect_equal(result$cov_theta, 1)
#    expect_equal(result$cov_psi, mm)
#    expect_equal(result$m_phi, 1)
#    expect_equal(result$m_theta, 2)
#    expect_equal(result$m_psi, 3:4)
#
#    # Discrete effect
#    result <- set_modargs(~ 1, ~ 1, ~ cov2, data)
#    mm     <- model.matrix(as.formula(~ cov2), data@site_cov)
#    expect_equal(result$M, 4)
#    expect_equal(result$code, "112000")
#    expect_equal(result$cov_phi, 1)
#    expect_equal(result$cov_theta, 1)
#    expect_equal(result$cov_psi, mm)
#    expect_equal(result$m_phi, 1)
#    expect_equal(result$m_theta, 2)
#    expect_equal(result$m_psi, 3:4)
#
#    # Additive effects
#    result <- set_modargs(~ 1, ~ 1, ~ cov1 + cov2, data)
#    mm     <- model.matrix(as.formula(~ cov1 + cov2), data@site_cov)
#    expect_equal(result$M, 5)
#    expect_equal(result$code, "112000")
#    expect_equal(result$cov_phi, 1)
#    expect_equal(result$cov_theta, 1)
#    expect_equal(result$cov_psi, mm)
#    expect_equal(result$m_phi, 1)
#    expect_equal(result$m_theta, 2)
#    expect_equal(result$m_psi, 3:5)
#
#    # Interactive effects
#    result <- set_modargs(~ 1, ~ 1, ~ cov1 * cov2, data)
#    mm     <- model.matrix(as.formula(~ cov1 * cov2), data@site_cov)
#    expect_equal(result$M, 6)
#    expect_equal(result$code, "112000")
#    expect_equal(result$cov_phi, 1)
#    expect_equal(result$cov_theta, 1)
#    expect_equal(result$cov_psi, mm)
#    expect_equal(result$m_phi, 1)
#    expect_equal(result$m_theta, 2)
#    expect_equal(result$m_psi, 3:6)
#})
#test_that("Setup for varying-theta/phi model works", {
#    cov1 <- rnorm(2); cov2 <- factor(1:2)
#    data <- occumbData(y = array(0, dim = rep(2, 3)),
#                       site_cov = list(cov1 = cov1,
#                                       cov2 = cov2))
#    expect_error(set_modargs(~ 1, ~ cov1, ~ 1, data),
#                 "Covariates in theta_formula are not yet supported")
#    expect_error(set_modargs(~ cov1, ~ cov1, ~ 1, data),
#                 "Covariates in phi_formula are not yet supported")
#})


### Set test cases -------------------------------------------------------------
phi <- theta <- c("i", "ij", "ijk")
psi <- c("i", "ij")
phi_shared <- theta_shared <- psi_shared <- c(FALSE, TRUE)

cases <- expand.grid(phi, theta, psi, phi_shared, theta_shared, psi_shared)
colnames(cases) <- c("phi", "theta", "psi",
                     "phi_shared", "theta_shared", "psi_shared")

### Tests for write_jags_model() -----------------------------------------------
test_that("JAGS code is correct for 144 available models", {
    for (i in 1:nrow(cases)) {
        ans <- readLines(system.file("jags",
                                     "occumb_template1.jags",
                                     package = "occumb"))

        if (cases$phi[i] == "i")
            ans <- c(ans,
                     "                x[i, j, k] ~ dgamma(phi[i], 1)")
        if (cases$phi[i] == "ij")
            ans <- c(ans,
                     "                x[i, j, k] ~ dgamma(phi[i, j], 1)")
        if (cases$phi[i] == "ijk")
            ans <- c(ans,
                     "                x[i, j, k] ~ dgamma(phi[i, j, k], 1)")

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template2.jags",
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

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template3.jags",
                                            package = "occumb")))

        if (cases$psi[i] == "i")
            ans <- c(ans,
                     "            z[i, j] ~ dbern(psi[i])")
        if (cases$psi[i] == "ij")
            ans <- c(ans,
                     "            z[i, j] ~ dbern(psi[i, j])")

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template4.jags",
                                            package = "occumb")))

        if (cases$phi_shared[i]) {
            if (cases$phi[i] == "i")
                ans <- c(ans,
                         "        log(phi[i]) <- inprod(alpha[i, ], cov_phi[]) + inprod(alpha_shared[], cov_phi_shared[i, ])")
            if (cases$phi[i] == "ij")
                ans <- c(ans,
                         "        log(phi[i, j]) <- inprod(alpha[i, ], cov_phi[i, j, ]) + inprod(alpha_shared[], cov_phi_shared[i, j, ])")
            if (cases$phi[i] == "ijk")
                ans <- c(ans,
                         "        log(phi[i, j, k]) <- inprod(alpha[i, ], cov_phi[i, j, k, ]) + inprod(alpha_shared[], cov_phi_shared[i, j, k, ])")
        } else {
            if (cases$phi[i] == "i")
                ans <- c(ans,
                         "        log(phi[i]) <- inprod(alpha[i, ], cov_phi[])")
            if (cases$phi[i] == "ij")
                ans <- c(ans,
                         "        log(phi[i, j]) <- inprod(alpha[i, ], cov_phi[i, j, ])")
            if (cases$phi[i] == "ijk")
                ans <- c(ans,
                         "        log(phi[i, j, k]) <- inprod(alpha[i, ], cov_phi[i, j, k, ])")
        }

        if (cases$theta_shared[i]) {
            if (cases$theta[i] == "i")
                ans <- c(ans,
                         "        logit(theta[i]) <- inprod(beta[i, ], cov_theta[]) + inprod(beta_shared[], cov_theta_shared[i, ])")
            if (cases$theta[i] == "ij")
                ans <- c(ans,
                         "        logit(theta[i, j]) <- inprod(beta[i, ], cov_theta[i, j, ]) + inprod(beta_shared[], cov_theta_shared[i, j, ])")
            if (cases$theta[i] == "ijk")
                ans <- c(ans,
                         "        logit(theta[i, j, k]) <- inprod(beta[i, ], cov_theta[i, j, k, ]) + inprod(beta_shared[], cov_theta_shared[i, j, k, ])")
        } else {
            if (cases$theta[i] == "i")
                ans <- c(ans,
                         "        logit(theta[i]) <- inprod(beta[i, ], cov_theta[])")
            if (cases$theta[i] == "ij")
                ans <- c(ans,
                         "        logit(theta[i, j]) <- inprod(beta[i, ], cov_theta[i, j, ])")
            if (cases$theta[i] == "ijk")
                ans <- c(ans,
                         "        logit(theta[i, j, k]) <- inprod(beta[i, ], cov_theta[i, j, k, ])")
        }

        if (cases$psi_shared[i]) {
            if (cases$psi[i] == "i")
                ans <- c(ans,
                         "        logit(psi[i]) <- inprod(gamma[i, ], cov_psi[]) + inprod(gamma_shared[], cov_psi_shared[i, ])")
            if (cases$psi[i] == "ij")
                ans <- c(ans,
                         "        logit(psi[i, j]) <- inprod(gamma[i, ], cov_psi[i, j, ]) + inprod(gamma_shared[], cov_psi_shared[i, j, ])")
        } else {
            if (cases$psi[i] == "i")
                ans <- c(ans,
                         "        logit(psi[i]) <- inprod(gamma[i, ], cov_psi[])")
            if (cases$psi[i] == "ij")
                ans <- c(ans,
                         "        logit(psi[i, j]) <- inprod(gamma[i, ], cov_psi[i, j, ])")
        }

        ans <- c(ans, readLines(system.file("jags",
                                            "occumb_template5.jags",
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

        res <- write_jags_model(phi          = cases$phi[i],
                                theta        = cases$theta[i],
                                psi          = cases$psi[i],
                                phi_shared   = cases$phi_shared[i],
                                theta_shared = cases$theta_shared[i],
                                psi_shared   = cases$psi_shared[i])
        expect_equal(res, ans)
    }
})

### Tests for set_data() -------------------------------------------------------
test_that("Data list is correct for 144 available models", {
    const <- list(I = sample.int(1E3, 1),
                  J = sample.int(1E3, 1),
                  K = sample.int(1E3, 1),
                  N = sample.int(1E3, 1),
                  y = sample.int(1E3, 1))
    prior_prec <- rnorm(1)
    prior_ulim <- rnorm(1)

    for (i in 1:nrow(cases)) {
        margs <- list(cov_phi          = rnorm(1),
                      cov_theta        = rnorm(1),
                      cov_psi          = rnorm(1),
                      cov_phi_shared   = rnorm(1),
                      cov_theta_shared = rnorm(1),
                      cov_psi_shared   = rnorm(1),
                      M                = sample.int(1E3, 1),
                      M_phi_shared     = sample.int(1E3, 1),
                      M_theta_shared   = sample.int(1E3, 1),
                      M_psi_shared     = sample.int(1E3, 1),
                      m_phi            = sample.int(1E3, 1),
                      m_theta          = sample.int(1E3, 1),
                      m_psi            = sample.int(1E3, 1),
                      phi_shared       = cases$phi_shared[i],
                      theta_shared     = cases$theta_shared[i],
                      psi_shared       = cases$psi_shared[i])

        ans <- list(I          = const$I,
                    J          = const$J,
                    K          = const$K,
                    N          = const$N,
                    y          = const$y,
                    cov_phi    = margs$cov_phi,
                    cov_theta  = margs$cov_theta,
                    cov_psi    = margs$cov_psi,
                    M          = margs$M,
                    m_phi      = margs$m_phi,
                    m_theta    = margs$m_theta,
                    m_psi      = margs$m_psi,
                    prior_prec = prior_prec,
                    prior_ulim = prior_ulim)

        if (cases$phi_shared[i])
            ans <- c(ans,
                     cov_phi_shared = margs$cov_phi_shared,
                     M_phi_shared   = margs$M_phi_shared)
        if (cases$theta_shared[i])
            ans <- c(ans,
                     cov_theta_shared = margs$cov_theta_shared,
                     M_theta_shared   = margs$M_theta_shared)
        if (cases$psi_shared[i])
            ans <- c(ans,
                     cov_psi_shared = margs$cov_psi_shared,
                     M_psi_shared   = margs$M_psi_shared)

        res <- set_data(const, margs, prior_prec, prior_ulim)
        expect_equal(res, ans)
    }
})

### Tests for inits_code() -----------------------------------------------------
test_that("Code for initial value function is correct for 144 available models", {
    const <- list(I = sample.int(1E3, 1),
                  J = sample.int(1E3, 1),
                  K = sample.int(1E3, 1),
                  N = sample.int(1E3, 1),
                  y = sample.int(1E3, 1))
    prior_prec <- rnorm(1)
    prior_ulim <- rnorm(1)

    for (i in 1:nrow(cases)) {
        margs <- list(cov_phi          = rnorm(1),
                      cov_theta        = rnorm(1),
                      cov_psi          = rnorm(1),
                      cov_phi_shared   = rnorm(1),
                      cov_theta_shared = rnorm(1),
                      cov_psi_shared   = rnorm(1),
                      M                = sample.int(1E3, 1),
                      M_phi_shared     = sample.int(1E3, 1),
                      M_theta_shared   = sample.int(1E3, 1),
                      M_psi_shared     = sample.int(1E3, 1),
                      m_phi            = sample.int(1E3, 1),
                      m_theta          = sample.int(1E3, 1),
                      m_psi            = sample.int(1E3, 1),
                      phi_shared       = cases$phi_shared[i],
                      theta_shared     = cases$theta_shared[i],
                      psi_shared       = cases$psi_shared[i])

        ans <- c(
            "function() {",
            "    list(z = matrix(1, const$I, const$J),",
            "         u = array(1, dim = c(const$I, const$J, const$K)),",
            "         x = array(stats::rnorm(const$I * const$J * const$K,",
            "                                mean = 1, sd = 0.1),",
            "                   dim = c(const$I, const$J, const$K)),")

        if (cases$phi_shared[i])
            ans <- c(ans, "         alpha_shared = stats::rnorm(margs$M_phi_shared, sd = 0.1),")
        if (cases$theta_shared[i])
            ans <- c(ans, "         beta_shared  = stats::rnorm(margs$M_theta_shared, sd = 0.1),")
        if (cases$psi_shared[i])
            ans <- c(ans, "         gamma_shared = stats::rnorm(margs$M_psi_shared, sd = 0.1),")

        ans <- c(ans,
            "         spec_eff = matrix(stats::rnorm(const$I * margs$M, sd = 0.1),",
            "                           const$I, margs$M),",
            "         Mu       = stats::rnorm(margs$M, sd = 0.1),",
            "         sigma    = stats::rnorm(margs$M, mean = 1, sd = 0.1),",
            "         rho      = matrix(stats::rnorm(margs$M^2, sd = 0.1),",
            "                           margs$M, margs$M))",
            "}")

        res <- inits_code(const, margs)
        expect_equal(res, ans)

    }
})

### Tests for set_params_monitored() -------------------------------------------
test_that("Parameter list is correct for 144 available models", {
    for (i in 1:nrow(cases)) {
        ans <- c("Mu", "sigma", "rho", "alpha", "beta", "gamma")

        if (cases$phi_shared[i])
            ans <- c(ans, "alpha_shared")
        if (cases$theta_shared[i])
            ans <- c(ans, "beta_shared")
        if (cases$psi_shared[i])
            ans <- c(ans, "gamma_shared")

        ans <- c(ans, c("phi", "theta", "psi", "z", "pi"))

        res <- set_params_monitored(phi_shared   = cases$phi_shared[i],
                                    theta_shared = cases$theta_shared[i],
                                    psi_shared   = cases$psi_shared[i])
        expect_equal(res, ans)
    }
})

