### Test data ------------------------------------------------------------------
I <- 2; J <- 3; K <- 4
cov1 <- rnorm(I)
cov2 <- factor(1:I)
cov3 <- rnorm(J)
cov4 <- factor(1:J)
cov5 <- matrix(rnorm(J * K), nrow = J)
cov6 <- matrix(rep(factor(1:K), each = J), nrow = J)
#cov2 <- as.character(1:I)
#cov4 <- as.character(factor(1:J))
#cov6 <- matrix(rep(as.character(1:K), each = J), nrow = J)
data <- occumbData(y = array(0, dim = c(I, J, K)),
                   spec_cov = list(cov1 = cov1, cov2 = cov2),
                   site_cov = list(cov3 = cov3, cov4 = cov4),
                   repl_cov = list(cov5 = cov5, cov6 = cov6))

### Tests for formula_phi_shared -----------------------------------------------
test_that("Temp: phi_shared correct", {
    ## Output
    # spec_cov (continuous)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ cov1, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "i")
    expect_true(result$phi_shared)
    expect_equal(result$M_phi_shared, 1)
    ans_cov <- array(dim = c(I, 1))
    for (i in 1:I) {
        ans_cov[i, 1] <- cov1[i]
    }
    colnames(ans_cov) <- "cov1"
    expect_equal(result$cov_phi_shared, ans_cov)

    # spec_cov (factor)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ cov2, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "i")
    expect_true(result$phi_shared)
    expect_equal(result$M_phi_shared, 1)
    ans_cov <- array(dim = c(I, 1))
    for (i in 1:I) {
        ans_cov[i, 1] <- as.numeric((cov2[i] == 2))
    }
    colnames(ans_cov) <- c("cov22")
    expect_equal(result$cov_phi_shared, ans_cov)

    # spec_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ cov1 * cov2, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "i")
    expect_true(result$phi_shared)
    expect_equal(result$M_phi_shared, 3)
    ans_cov <- array(dim = c(I, 3))
    for (i in 1:I) {
        ans_cov[i, 1] <- cov1[i]
        ans_cov[i, 2] <- as.numeric((cov2[i] == 2))
        ans_cov[i, 3] <- cov1[i] * as.numeric((cov2[i] == 2))
    }
    colnames(ans_cov) <- c("cov1", "cov22", "cov1:cov22")
    expect_equal(result$cov_phi_shared, ans_cov)

    # site_cov (continuous)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ cov3, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "ij")
    expect_true(result$phi_shared)
    expect_equal(result$M_phi_shared, 1)
    ans_cov <- array(dim = c(I, J, 1))
    for (i in 1:I) {
        for (j in 1:J) {
            ans_cov[i, j, 1] <- cov3[j]
        }
    }
    dimnames(ans_cov)[[3]] <- "cov3"
    expect_equal(result$cov_phi_shared, ans_cov)

    # site_cov (factor)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ cov4, ~ 1, ~ 1, data = data)
    expect_true(result$phi_shared)
    expect_equal(result$M_phi_shared, 2)
    ans_cov <- array(dim = c(I, J, 2))
    for (i in 1:I) {
        for (j in 1:J) {
            ans_cov[i, j, 1] <- as.numeric(cov4[j] == 2)
            ans_cov[i, j, 2] <- as.numeric(cov4[j] == 3)
        }
    }
    dimnames(ans_cov)[[3]] <- c("cov42", "cov43")
    expect_equal(result$cov_phi_shared, ans_cov)

    # site_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ cov3 * cov4, ~ 1, ~ 1, data = data)
    expect_true(result$phi_shared)
    expect_equal(result$M_phi_shared, 5)
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
    expect_equal(result$cov_phi_shared, ans_cov)

    # spec_cov * site_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ cov1 * cov4, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "ij")
    expect_true(result$phi_shared)
    expect_equal(result$M_phi_shared, 5)
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
    expect_equal(result$cov_phi_shared, ans_cov)

    # repl_cov (continuous)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ cov5, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "ijk")
    expect_true(result$phi_shared)
    expect_equal(result$M_phi_shared, 1)
    ans_cov <- array(dim = c(I, J, K, 1))
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                ans_cov[i, j, k, 1] <- cov5[j, k]
            }
        }
    }
    dimnames(ans_cov)[[4]] <- "cov5"
    expect_equal(result$cov_phi_shared, ans_cov)

    # repl_cov (factor)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ cov6, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "ijk")
    expect_true(result$phi_shared)
    expect_equal(result$M_phi_shared, 3)
    ans_cov <- array(dim = c(I, J, K, 3))
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                ans_cov[i, j, k, 1] <- as.numeric(cov6[j, k] == 2)
                ans_cov[i, j, k, 2] <- as.numeric(cov6[j, k] == 3)
                ans_cov[i, j, k, 3] <- as.numeric(cov6[j, k] == 4)
            }
        }
    }
    dimnames(ans_cov)[[4]] <- c("cov62", "cov63", "cov64")
    expect_equal(result$cov_phi_shared, ans_cov)

    # repl_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ cov5 * cov6, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "ijk")
    expect_true(result$phi_shared)
    expect_equal(result$M_phi_shared, 7)
    ans_cov <- array(dim = c(I, J, K, 7))
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                ans_cov[i, j, k, 1] <- cov5[j, k]
                ans_cov[i, j, k, 2] <- as.numeric(cov6[j, k] == 2)
                ans_cov[i, j, k, 3] <- as.numeric(cov6[j, k] == 3)
                ans_cov[i, j, k, 4] <- as.numeric(cov6[j, k] == 4)
                ans_cov[i, j, k, 5] <- cov5[j, k] * as.numeric(cov6[j, k] == 2)
                ans_cov[i, j, k, 6] <- cov5[j, k] * as.numeric(cov6[j, k] == 3)
                ans_cov[i, j, k, 7] <- cov5[j, k] * as.numeric(cov6[j, k] == 4)
            }
        }
    }
    dimnames(ans_cov)[[4]] <- c("cov5", "cov62", "cov63", "cov64",
                                "cov5:cov62", "cov5:cov63", "cov5:cov64")
    expect_equal(result$cov_phi_shared, ans_cov)

    # spec_cov * site_cov * repl_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ cov2 * cov3 * cov5, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "ijk")
    expect_true(result$phi_shared)
    expect_equal(result$M_phi_shared, 7)
    ans_cov <- array(dim = c(I, J, K, 7))
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                ans_cov[i, j, k, 1] <- as.numeric(cov2[i] == 2)
                ans_cov[i, j, k, 2] <- cov3[j]
                ans_cov[i, j, k, 3] <- cov5[j, k]
                ans_cov[i, j, k, 4] <- as.numeric(cov2[i] == 2) * cov3[j]
                ans_cov[i, j, k, 5] <- as.numeric(cov2[i] == 2) * cov5[j, k]
                ans_cov[i, j, k, 6] <- cov3[j] * cov5[j, k]
                ans_cov[i, j, k, 7] <- as.numeric(cov2[i] == 2) * cov3[j] * cov5[j, k]
            }
        }
    }
    dimnames(ans_cov)[[4]] <- c("cov22", "cov3", "cov5",
                                "cov22:cov3", "cov22:cov5", "cov3:cov5",
                                "cov22:cov3:cov5")
    expect_equal(result$cov_phi_shared, ans_cov)

    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ cov1 * cov4 * cov5, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "ijk")
    expect_true(result$phi_shared)
    expect_equal(result$M_phi_shared, 11)
    ans_cov <- array(dim = c(I, J, K, 11))
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                ans_cov[i, j, k, 1] <- cov1[i]
                ans_cov[i, j, k, 2] <- as.numeric(cov4[j] == 2)
                ans_cov[i, j, k, 3] <- as.numeric(cov4[j] == 3)
                ans_cov[i, j, k, 4] <- cov5[j, k]
                ans_cov[i, j, k, 5] <- cov1[i] * as.numeric(cov4[j] == 2)
                ans_cov[i, j, k, 6] <- cov1[i] * as.numeric(cov4[j] == 3)
                ans_cov[i, j, k, 7] <- cov1[i] * cov5[j, k]
                ans_cov[i, j, k, 8] <- as.numeric(cov4[j] == 2) * cov5[j, k]
                ans_cov[i, j, k, 9] <- as.numeric(cov4[j] == 3) * cov5[j, k]
                ans_cov[i, j, k, 10] <- cov1[i] * as.numeric(cov4[j] == 2) * cov5[j, k]
                ans_cov[i, j, k, 11] <- cov1[i] * as.numeric(cov4[j] == 3) * cov5[j, k]
            }
        }
    }
    dimnames(ans_cov)[[4]] <- c("cov1", "cov42", "cov43", "cov5",
                                "cov1:cov42", "cov1:cov43",
                                "cov1:cov5",
                                "cov42:cov5", "cov43:cov5",
                                "cov1:cov42:cov5", "cov1:cov43:cov5")
    expect_equal(result$cov_phi_shared, ans_cov)

    ## Errors and Warnings
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, ~ 0, ~ 1, ~ 1, data = data),
                 "No intercept in formula_phi_shared: remove 0 or -1 from the formula")
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, ~ -1, ~ 1, ~ 1, data = data),
                 "No intercept in formula_phi_shared: remove 0 or -1 from the formula")
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, ~ xxx, ~ 1, ~ 1, data = data),
                             sprintf("Unexpected terms in formula_phi_shared: %s
Make sure they are found in either spec_cov, site_cov, or repl_cov.", "xxx"))
})

### Tests for formula_theta_shared ---------------------------------------------
test_that("Temp: theta_shared correct", {
    ## Output
    # spec_cov (continuous)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ cov1, ~ 1, data = data)
    expect_equal(result$theta, "i")
    expect_true(result$theta_shared)
    expect_equal(result$M_theta_shared, 1)
    ans_cov <- array(dim = c(I, 1))
    for (i in 1:I) {
        ans_cov[i, 1] <- cov1[i]
    }
    colnames(ans_cov) <- "cov1"
    expect_equal(result$cov_theta_shared, ans_cov)

    # spec_cov (factor)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ cov2, ~ 1, data = data)
    expect_equal(result$theta, "i")
    expect_true(result$theta_shared)
    expect_equal(result$M_theta_shared, 1)
    ans_cov <- array(dim = c(I, 1))
    for (i in 1:I) {
        ans_cov[i, 1] <- as.numeric((cov2[i] == 2))
    }
    colnames(ans_cov) <- c("cov22")
    expect_equal(result$cov_theta_shared, ans_cov)

    # spec_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ cov1 * cov2, ~ 1, data = data)
    expect_equal(result$theta, "i")
    expect_true(result$theta_shared)
    expect_equal(result$M_theta_shared, 3)
    ans_cov <- array(dim = c(I, 3))
    for (i in 1:I) {
        ans_cov[i, 1] <- cov1[i]
        ans_cov[i, 2] <- as.numeric((cov2[i] == 2))
        ans_cov[i, 3] <- cov1[i] * as.numeric((cov2[i] == 2))
    }
    colnames(ans_cov) <- c("cov1", "cov22", "cov1:cov22")
    expect_equal(result$cov_theta_shared, ans_cov)

    # site_cov (continuous)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ cov3, ~ 1, data = data)
    expect_equal(result$theta, "ij")
    expect_true(result$theta_shared)
    expect_equal(result$M_theta_shared, 1)
    ans_cov <- array(dim = c(I, J, 1))
    for (i in 1:I) {
        for (j in 1:J) {
            ans_cov[i, j, 1] <- cov3[j]
        }
    }
    dimnames(ans_cov)[[3]] <- "cov3"
    expect_equal(result$cov_theta_shared, ans_cov)

    # site_cov (factor)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ cov4, ~ 1, data = data)
    expect_true(result$theta_shared)
    expect_equal(result$M_theta_shared, 2)
    ans_cov <- array(dim = c(I, J, 2))
    for (i in 1:I) {
        for (j in 1:J) {
            ans_cov[i, j, 1] <- as.numeric(cov4[j] == 2)
            ans_cov[i, j, 2] <- as.numeric(cov4[j] == 3)
        }
    }
    dimnames(ans_cov)[[3]] <- c("cov42", "cov43")
    expect_equal(result$cov_theta_shared, ans_cov)

    # site_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ cov3 * cov4, ~ 1, data = data)
    expect_true(result$theta_shared)
    expect_equal(result$M_theta_shared, 5)
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
    expect_equal(result$cov_theta_shared, ans_cov)

    # spec_cov * site_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ cov1 * cov4, ~ 1, data = data)
    expect_equal(result$theta, "ij")
    expect_true(result$theta_shared)
    expect_equal(result$M_theta_shared, 5)
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
    expect_equal(result$cov_theta_shared, ans_cov)

    # repl_cov (continuous)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ cov5, ~ 1, data = data)
    expect_equal(result$theta, "ijk")
    expect_true(result$theta_shared)
    expect_equal(result$M_theta_shared, 1)
    ans_cov <- array(dim = c(I, J, K, 1))
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                ans_cov[i, j, k, 1] <- cov5[j, k]
            }
        }
    }
    dimnames(ans_cov)[[4]] <- "cov5"
    expect_equal(result$cov_theta_shared, ans_cov)

    # repl_cov (factor)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ cov6, ~ 1, data = data)
    expect_equal(result$theta, "ijk")
    expect_true(result$theta_shared)
    expect_equal(result$M_theta_shared, 3)
    ans_cov <- array(dim = c(I, J, K, 3))
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                ans_cov[i, j, k, 1] <- as.numeric(cov6[j, k] == 2)
                ans_cov[i, j, k, 2] <- as.numeric(cov6[j, k] == 3)
                ans_cov[i, j, k, 3] <- as.numeric(cov6[j, k] == 4)
            }
        }
    }
    dimnames(ans_cov)[[4]] <- c("cov62", "cov63", "cov64")
    expect_equal(result$cov_theta_shared, ans_cov)

    # repl_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ cov5 * cov6, ~ 1, data = data)
    expect_equal(result$theta, "ijk")
    expect_true(result$theta_shared)
    expect_equal(result$M_theta_shared, 7)
    ans_cov <- array(dim = c(I, J, K, 7))
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                ans_cov[i, j, k, 1] <- cov5[j, k]
                ans_cov[i, j, k, 2] <- as.numeric(cov6[j, k] == 2)
                ans_cov[i, j, k, 3] <- as.numeric(cov6[j, k] == 3)
                ans_cov[i, j, k, 4] <- as.numeric(cov6[j, k] == 4)
                ans_cov[i, j, k, 5] <- cov5[j, k] * as.numeric(cov6[j, k] == 2)
                ans_cov[i, j, k, 6] <- cov5[j, k] * as.numeric(cov6[j, k] == 3)
                ans_cov[i, j, k, 7] <- cov5[j, k] * as.numeric(cov6[j, k] == 4)
            }
        }
    }
    dimnames(ans_cov)[[4]] <- c("cov5", "cov62", "cov63", "cov64",
                                "cov5:cov62", "cov5:cov63", "cov5:cov64")
    expect_equal(result$cov_theta_shared, ans_cov)

    # spec_cov * site_cov * repl_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ cov2 * cov3 * cov5, ~ 1, data = data)
    expect_equal(result$theta, "ijk")
    expect_true(result$theta_shared)
    expect_equal(result$M_theta_shared, 7)
    ans_cov <- array(dim = c(I, J, K, 7))
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                ans_cov[i, j, k, 1] <- as.numeric(cov2[i] == 2)
                ans_cov[i, j, k, 2] <- cov3[j]
                ans_cov[i, j, k, 3] <- cov5[j, k]
                ans_cov[i, j, k, 4] <- as.numeric(cov2[i] == 2) * cov3[j]
                ans_cov[i, j, k, 5] <- as.numeric(cov2[i] == 2) * cov5[j, k]
                ans_cov[i, j, k, 6] <- cov3[j] * cov5[j, k]
                ans_cov[i, j, k, 7] <- as.numeric(cov2[i] == 2) * cov3[j] * cov5[j, k]
            }
        }
    }
    dimnames(ans_cov)[[4]] <- c("cov22", "cov3", "cov5",
                                "cov22:cov3", "cov22:cov5", "cov3:cov5",
                                "cov22:cov3:cov5")
    expect_equal(result$cov_theta_shared, ans_cov)

    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ cov1 * cov4 * cov5, ~ 1, data = data)
    expect_equal(result$theta, "ijk")
    expect_true(result$theta_shared)
    expect_equal(result$M_theta_shared, 11)
    ans_cov <- array(dim = c(I, J, K, 11))
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                ans_cov[i, j, k, 1] <- cov1[i]
                ans_cov[i, j, k, 2] <- as.numeric(cov4[j] == 2)
                ans_cov[i, j, k, 3] <- as.numeric(cov4[j] == 3)
                ans_cov[i, j, k, 4] <- cov5[j, k]
                ans_cov[i, j, k, 5] <- cov1[i] * as.numeric(cov4[j] == 2)
                ans_cov[i, j, k, 6] <- cov1[i] * as.numeric(cov4[j] == 3)
                ans_cov[i, j, k, 7] <- cov1[i] * cov5[j, k]
                ans_cov[i, j, k, 8] <- as.numeric(cov4[j] == 2) * cov5[j, k]
                ans_cov[i, j, k, 9] <- as.numeric(cov4[j] == 3) * cov5[j, k]
                ans_cov[i, j, k, 10] <- cov1[i] * as.numeric(cov4[j] == 2) * cov5[j, k]
                ans_cov[i, j, k, 11] <- cov1[i] * as.numeric(cov4[j] == 3) * cov5[j, k]
            }
        }
    }
    dimnames(ans_cov)[[4]] <- c("cov1", "cov42", "cov43", "cov5",
                                "cov1:cov42", "cov1:cov43",
                                "cov1:cov5",
                                "cov42:cov5", "cov43:cov5",
                                "cov1:cov42:cov5", "cov1:cov43:cov5")
    expect_equal(result$cov_theta_shared, ans_cov)

    ## Errors and Warnings
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 0, ~ 1, data = data),
                 "No intercept in formula_theta_shared: remove 0 or -1 from the formula")
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ -1, ~ 1, data = data),
                 "No intercept in formula_theta_shared: remove 0 or -1 from the formula")
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ xxx, ~ 1, data = data),
                             sprintf("Unexpected terms in formula_theta_shared: %s
Make sure they are found in either spec_cov, site_cov, or repl_cov.", "xxx"))
})

### Tests for formula_psi_shared -----------------------------------------------
test_that("Temp: psi_shared correct", {
    ## Output
    # spec_cov (continuous)
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ cov1, data = data)
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
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ cov2, data = data)
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
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ cov1 * cov2, data = data)
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
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ cov3, data = data)
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
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ cov4, data = data)
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
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ cov3 * cov4, data = data)
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
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ cov1 * cov4, data = data)
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
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ cov5, data = data),
                             sprintf("Unexpected terms in formula_psi_shared: %s
Note that only site covariates, species covariates, or their interactions are allowed for formula_psi_shared.", "cov5"))

    ## Errors and Warnings
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ 0, data = data),
                 sprintf("No intercept in formula_%s: remove 0 or -1 from the formula", "psi_shared"))
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ -1, data = data),
                 sprintf("No intercept in formula_%s: remove 0 or -1 from the formula", "psi_shared"))
    expect_error(set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ xxx, data = data),
                             sprintf("Unexpected terms in formula_psi_shared: %s
Note that only site covariates, species covariates, or their interactions are allowed for formula_psi_shared.", "xxx"))
})

### Tests for formula_phi ------------------------------------------------------
test_that("Arguments are correct for the null model", {
    ## Output
    # null model
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "i")
    expect_equal(result$M, 3)
    expect_equal(result$cov_phi, 1)
    expect_equal(result$m_phi, 1)
})

test_that("Errors are correct for the spec_cov in phi", {
    # spec_cov (not allowed)
    expect_error(set_modargs(~ cov1, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, data = data),
                 sprintf("Unexpected terms in formula_phi: %s
Note that species covariates are not allowed for formula_phi.",
                         "cov1"))
})

test_that("Temp: phi correct", {
    # site_cov (continuous)
    result <- set_modargs(~ cov3, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "ij")
    expect_equal(result$M, 4)
    ans_cov <- array(dim = c(I, J, 2))
    for (i in 1:I) {
        for (j in 1:J) {
            ans_cov[i, j, 1] <- 1
            ans_cov[i, j, 2] <- cov3[j]
        }
    }
    dimnames(ans_cov)[[3]] <- c("(Intercept)", "cov3")
    expect_equal(result$cov_phi, ans_cov)
    expect_equal(result$m_phi, 1:2)

    # site_cov (factor)
    result <- set_modargs(~ cov4, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "ij")
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
    expect_equal(result$cov_phi, ans_cov)
    expect_equal(result$m_phi, 1:3)

    # site_cov (interaction)
    result <- set_modargs(~ cov3 * cov4, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "ij")
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
    expect_equal(result$cov_phi, ans_cov)
    expect_equal(result$m_phi, 1:6)

    # repl_cov (continuous)
    result <- set_modargs(~ cov5, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "ijk")
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
    expect_equal(result$cov_phi, ans_cov)
    expect_equal(result$m_phi, 1:2)

    # repl_cov (factor)
    result <- set_modargs(~ cov6, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "ijk")
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
    expect_equal(result$cov_phi, ans_cov)
    expect_equal(result$m_phi, 1:4)

    # repl_cov (interaction)
    result <- set_modargs(~ cov5 * cov6, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "ijk")
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
    expect_equal(result$cov_phi, ans_cov)
    expect_equal(result$m_phi, 1:8)

    # site_cov * repl_cov (interaction)
    result <- set_modargs(~ cov3 * cov6, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$phi, "ijk")
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
    expect_equal(result$cov_phi, ans_cov)
    expect_equal(result$m_phi, 1:8)

    ## Errors and Warnings
    expect_error(set_modargs(~ 0, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, data = data),
                 sprintf("No intercept in formula_%s: remove 0 or -1 from the formula", "phi"))
    expect_error(set_modargs(~ -1, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, data = data),
                 sprintf("No intercept in formula_%s: remove 0 or -1 from the formula", "phi"))
    expect_error(set_modargs(~ xxx, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, data = data),
                 sprintf("Unexpected terms in formula_phi: %s
Note that species covariates are not allowed for formula_phi",
                         "xxx"))
})

### Tests for formula_theta ----------------------------------------------------
test_that("Temp: theta correct", {
    ## Output
    # null model
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$theta, "i")
    expect_equal(result$M, 3)
    expect_equal(result$cov_theta, 1)
    expect_equal(result$m_theta, 2)

    # spec_cov (not allowed)
    expect_error(set_modargs(~ 1, ~ cov1, ~ 1, ~ 1, ~ 1, ~ 1, data = data),
                 sprintf("Unexpected terms in formula_theta: %s
Note that species covariates are not allowed for formula_theta.",
                         "cov1"))

    # site_cov (continuous)
    result <- set_modargs(~ 1, ~ cov3, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$theta, "ij")
    expect_equal(result$M, 4)
    ans_cov <- array(dim = c(J, 2))
    for (j in 1:J) {
        ans_cov[j, 1] <- 1
        ans_cov[j, 2] <- cov3[j]
    }
    colnames(ans_cov) <- c("(Intercept)", "cov3")
    expect_equal(result$cov_theta, ans_cov)
    expect_equal(result$m_theta, 2:3)

    # site_cov (factor)
    result <- set_modargs(~ 1, ~ cov4, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$theta, "ij")
    expect_equal(result$M, 5)
    ans_cov <- array(dim = c(J, 3))
    for (j in 1:J) {
        ans_cov[j, 1] <- 1
        ans_cov[j, 2] <- as.numeric((cov4[j] == 2))
        ans_cov[j, 3] <- as.numeric((cov4[j] == 3))
    }
    colnames(ans_cov) <- c("(Intercept)", "cov42", "cov43")
    expect_equal(result$cov_theta, ans_cov)
    expect_equal(result$m_theta, 2:4)

    # site_cov (interaction)
    result <- set_modargs(~ 1, ~ cov3 * cov4, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$theta, "ij")
    expect_equal(result$M, 8)
    ans_cov <- array(dim = c(J, 6))
    for (j in 1:J) {
        ans_cov[j, 1] <- 1
        ans_cov[j, 2] <- cov3[j]
        ans_cov[j, 3] <- as.numeric(cov4[j] == 2)
        ans_cov[j, 4] <- as.numeric(cov4[j] == 3)
        ans_cov[j, 5] <- cov3[j] * as.numeric(cov4[j] == 2)
        ans_cov[j, 6] <- cov3[j] * as.numeric(cov4[j] == 3)
    }
    colnames(ans_cov) <- c("(Intercept)", "cov3", "cov42", "cov43",
                           "cov3:cov42", "cov3:cov43")
    expect_equal(result$cov_theta, ans_cov)
    expect_equal(result$m_theta, 2:7)

    # repl_cov (continuous)
    result <- set_modargs(~ 1, ~ cov5, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$theta, "ijk")
    expect_equal(result$M, 4)
    ans_cov <- array(dim = c(J, K, 2))
    for (j in 1:J) {
        for (k in 1:K) {
            ans_cov[j, k, 1] <- 1
            ans_cov[j, k, 2] <- cov5[j, k]
        }
    }
    dimnames(ans_cov)[[3]] <- c("(Intercept)", "cov5")
    expect_equal(result$cov_theta, ans_cov)
    expect_equal(result$m_theta, 2:3)

    # repl_cov (factor)
    result <- set_modargs(~ 1, ~ cov6, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$theta, "ijk")
    expect_equal(result$M, 6)
    ans_cov <- array(dim = c(J, K, 4))
    for (j in 1:J) {
        for (k in 1:K) {
            ans_cov[j, k, 1] <- 1
            ans_cov[j, k, 2] <- as.numeric((cov6[j, k] == 2))
            ans_cov[j, k, 3] <- as.numeric((cov6[j, k] == 3))
            ans_cov[j, k, 4] <- as.numeric((cov6[j, k] == 4))
        }
    }
    dimnames(ans_cov)[[3]] <- c("(Intercept)", "cov62", "cov63", "cov64")
    expect_equal(result$cov_theta, ans_cov)
    expect_equal(result$m_theta, 2:5)

    # repl_cov (interaction)
    result <- set_modargs(~ 1, ~ cov5 * cov6, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$theta, "ijk")
    expect_equal(result$M, 10)
    ans_cov <- array(dim = c(J, K, 8))
    for (j in 1:J) {
        for (k in 1:K) {
            ans_cov[j, k, 1] <- 1
            ans_cov[j, k, 2] <- cov5[j, k]
            ans_cov[j, k, 3] <- as.numeric(cov6[j, k] == 2)
            ans_cov[j, k, 4] <- as.numeric(cov6[j, k] == 3)
            ans_cov[j, k, 5] <- as.numeric(cov6[j, k] == 4)
            ans_cov[j, k, 6] <- cov5[j, k] * as.numeric(cov6[j, k] == 2)
            ans_cov[j, k, 7] <- cov5[j, k] * as.numeric(cov6[j, k] == 3)
            ans_cov[j, k, 8] <- cov5[j, k] * as.numeric(cov6[j, k] == 4)
        }
    }
    dimnames(ans_cov)[[3]] <- c("(Intercept)", "cov5", "cov62", "cov63", "cov64",
                                "cov5:cov62", "cov5:cov63", "cov5:cov64")
    expect_equal(result$cov_theta, ans_cov)
    expect_equal(result$m_theta, 2:9)

    # site_cov * repl_cov (interaction)
    result <- set_modargs(~ 1, ~ cov3 * cov6, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$theta, "ijk")
    expect_equal(result$M, 10)
    ans_cov <- array(dim = c(J, K, 8))
    for (j in 1:J) {
        for (k in 1:K) {
            ans_cov[j, k, 1] <- 1
            ans_cov[j, k, 2] <- cov3[j]
            ans_cov[j, k, 3] <- as.numeric(cov6[j, k] == 2)
            ans_cov[j, k, 4] <- as.numeric(cov6[j, k] == 3)
            ans_cov[j, k, 5] <- as.numeric(cov6[j, k] == 4)
            ans_cov[j, k, 6] <- cov3[j] * as.numeric(cov6[j, k] == 2)
            ans_cov[j, k, 7] <- cov3[j] * as.numeric(cov6[j, k] == 3)
            ans_cov[j, k, 8] <- cov3[j] * as.numeric(cov6[j, k] == 4)
        }
    }
    dimnames(ans_cov)[[3]] <- c("(Intercept)", "cov3", "cov62", "cov63", "cov64",
                                "cov3:cov62", "cov3:cov63", "cov3:cov64")
    expect_equal(result$cov_theta, ans_cov)
    expect_equal(result$m_theta, 2:9)

    ## Errors and Warnings
    expect_error(set_modargs(~ 1, ~ 0, ~ 1, ~ 1, ~ 1, ~ 1, data = data),
                 sprintf("No intercept in formula_%s: remove 0 or -1 from the formula", "theta"))
    expect_error(set_modargs(~ 1, ~ -1, ~ 1, ~ 1, ~ 1, ~ 1, data = data),
                 sprintf("No intercept in formula_%s: remove 0 or -1 from the formula", "theta"))
    expect_error(set_modargs(~ 1, ~ xxx, ~ 1, ~ 1, ~ 1, ~ 1, data = data),
                 sprintf("Unexpected terms in formula_theta: %s
Note that species covariates are not allowed for formula_theta.",
                         "xxx"))
})

### Tests for formula_psi ------------------------------------------------------
test_that("Temp: psi correct", {
    ## Output
    # null model
    result <- set_modargs(~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$psi, "i")
    expect_equal(result$M, 3)
    expect_equal(result$cov_psi, 1)
    expect_equal(result$m_psi, 3)

    # spec_cov (not allowed)
    expect_error(set_modargs(~ 1, ~ 1, ~ cov1, ~ 1, ~ 1, ~ 1, data = data),
                 sprintf("Unexpected terms in formula_psi: %s
Note that only site covariates are allowed for formula_psi.",
                         "cov1"))

    # site_cov (continuous)
    result <- set_modargs(~ 1, ~ 1, ~ cov3, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$psi, "ij")
    expect_equal(result$M, 4)
    ans_cov <- array(dim = c(J, 2))
    for (j in 1:J) {
        ans_cov[j, 1] <- 1
        ans_cov[j, 2] <- cov3[j]
    }
    colnames(ans_cov) <- c("(Intercept)", "cov3")
    expect_equal(result$cov_psi, ans_cov)
    expect_equal(result$m_psi, 3:4)

    # site_cov (factor)
    result <- set_modargs(~ 1, ~ 1, ~ cov4, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$psi, "ij")
    expect_equal(result$M, 5)
    ans_cov <- array(dim = c(J, 3))
    for (j in 1:J) {
        ans_cov[j, 1] <- 1
        ans_cov[j, 2] <- as.numeric(cov4[j] == 2)
        ans_cov[j, 3] <- as.numeric(cov4[j] == 3)
    }
    colnames(ans_cov) <- c("(Intercept)", "cov42", "cov43")
    expect_equal(result$cov_psi, ans_cov)
    expect_equal(result$m_psi, 3:5)

    # site_cov (interaction)
    result <- set_modargs(~ 1, ~ 1, ~ cov3 * cov4, ~ 1, ~ 1, ~ 1, data = data)
    expect_equal(result$psi, "ij")
    expect_equal(result$M, 8)
    ans_cov <- array(dim = c(J, 6))
    for (j in 1:J) {
        ans_cov[j, 1] <- 1
        ans_cov[j, 2] <- cov3[j]
        ans_cov[j, 3] <- as.numeric(cov4[j] == 2)
        ans_cov[j, 4] <- as.numeric(cov4[j] == 3)
        ans_cov[j, 5] <- cov3[j] * as.numeric(cov4[j] == 2)
        ans_cov[j, 6] <- cov3[j] * as.numeric(cov4[j] == 3)
    }
    colnames(ans_cov) <- c("(Intercept)", "cov3", "cov42", "cov43",
                           "cov3:cov42", "cov3:cov43")
    expect_equal(result$cov_psi, ans_cov)
    expect_equal(result$m_psi, 3:8)

    # repl_cov (not allowed)
    expect_error(set_modargs(~ 1, ~ 1, ~ cov5, ~ 1, ~ 1, ~ 1, data = data),
                 sprintf("Unexpected terms in formula_psi: %s
Note that only site covariates are allowed for formula_psi.",
                         "cov5"))

    ## Errors and Warnings
    expect_error(set_modargs(~ 1, ~ 1, ~ 0, ~ 1, ~ 1, ~ 1, data = data),
                 sprintf("No intercept in formula_%s: remove 0 or -1 from the formula", "psi"))
    expect_error(set_modargs(~ 1, ~ 1, ~ -1, ~ 1, ~ 1, ~ 1, data = data),
                 sprintf("No intercept in formula_%s: remove 0 or -1 from the formula", "psi"))
    expect_error(set_modargs(~ 1, ~ 1, ~ xxx, ~ 1, ~ 1, ~ 1, data = data),
                 sprintf("Unexpected terms in formula_psi: %s
Note that only site covariates are allowed for formula_psi.",
                         "xxx"))
})

