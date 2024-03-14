### Tests for eval_util_L/R ----------------------------------------------------
I <- 20 # Number of species
J <- 2 # Number of sites
K <- 2 # Number of replicates
data <- occumbData(
    y = array(sample.int(I * J * K), dim = c(I, J, K)),
    spec_cov = list(cov1 = rnorm(I)),
    site_cov = list(cov2 = rnorm(J),
                    cov3 = factor(1:J)),
    repl_cov = list(cov4 = matrix(rnorm(J * K), J, K)))

# Fitting a null model (includes only species-specific intercepts)
res0 <- occumb(data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)

# Add species-specific effects of site covariates in occupancy probabilities
res1 <- occumb(formula_psi = ~ cov2, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)
res1a <- occumb(formula_theta = ~ cov2, data = data,
                n.chains = 1, n.adapt = 0, n.burnin = 0,
                n.thin = 1, n.iter = 10, verbose = FALSE)
res1b <- occumb(formula_phi = ~ cov2, data = data,
                n.chains = 1, n.adapt = 0, n.burnin = 0,
                n.thin = 1, n.iter = 10, verbose = FALSE)
res2 <- occumb(formula_psi = ~ cov3, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)
res3 <- occumb(formula_psi = ~ cov2 * cov3, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)

# Add species covariate in the three parameters
res4 <- occumb(formula_phi_shared = ~ cov1, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)
res5 <- occumb(formula_theta_shared = ~ cov1, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)
res6 <- occumb(formula_psi_shared = ~ cov1, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)

# Add replicate covariates
res7 <- occumb(formula_phi = ~ cov4, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)
res8 <- occumb(formula_theta = ~ cov4, data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)

test_that("eval_util_L() outputs a data frame with the additional Utility column", {
    settings <- data.frame(K = rep(1, 3), N = rep(1, 3), x = NA)

    # Null model
    test0 <- eval_util_L(settings, res0, cores = 1)
    checkmate::expect_data_frame(test0)
    expect_equal(colnames(test0), c(colnames(settings), "Utility"))
    expect_equal(test0[, -ncol(test0)], settings)

    # Model with site covariates
    test1 <- eval_util_L(settings, res1, cores = 1)
    checkmate::expect_data_frame(test1)
    expect_equal(colnames(test1), c(colnames(settings), "Utility"))
    expect_equal(test1[, -ncol(test1)], settings)

    test1a <- eval_util_L(settings, res1a, cores = 1)
    checkmate::expect_data_frame(test1a)
    expect_equal(colnames(test1a), c(colnames(settings), "Utility"))
    expect_equal(test1a[, -ncol(test1a)], settings)

    test1b <- eval_util_L(settings, res1b, cores = 1)
    checkmate::expect_data_frame(test1b)
    expect_equal(colnames(test1b), c(colnames(settings), "Utility"))
    expect_equal(test1b[, -ncol(test1b)], settings)

    test2 <- eval_util_L(settings, res2, cores = 1)
    checkmate::expect_data_frame(test2)
    expect_equal(colnames(test2), c(colnames(settings), "Utility"))
    expect_equal(test2[, -ncol(test2)], settings)

    test3 <- eval_util_L(settings, res3, cores = 1)
    checkmate::expect_data_frame(test3)
    expect_equal(colnames(test3), c(colnames(settings), "Utility"))
    expect_equal(test3[, -ncol(test3)], settings)

    # Model with species covariates
    test4 <- eval_util_L(settings, res4, cores = 1)
    checkmate::expect_data_frame(test4)
    expect_equal(colnames(test4), c(colnames(settings), "Utility"))
    expect_equal(test4[, -ncol(test4)], settings)

    test5 <- eval_util_L(settings, res5, cores = 1)
    checkmate::expect_data_frame(test5)
    expect_equal(colnames(test5), c(colnames(settings), "Utility"))
    expect_equal(test5[, -ncol(test5)], settings)

    test6 <- eval_util_L(settings, res6, cores = 1)
    checkmate::expect_data_frame(test6)
    expect_equal(colnames(test6), c(colnames(settings), "Utility"))
    expect_equal(test6[, -ncol(test6)], settings)

    # Use z, theta, phi arguments
    test0_z <- eval_util_L(settings, res0,
                           z = array(1, c(1, I, J)),
                           cores = 1)
    checkmate::expect_data_frame(test0_z)
    expect_equal(colnames(test0_z), c(colnames(settings), "Utility"))
    expect_equal(test0_z[, -ncol(test0_z)], settings)

    test0_theta1 <- eval_util_L(settings, res0,
                                theta = matrix(runif(1 * I), 1, I),
                                cores = 1)
    checkmate::expect_data_frame(test0_theta1)
    expect_equal(colnames(test0_theta1), c(colnames(settings), "Utility"))
    expect_equal(test0_theta1[, -ncol(test0_theta1)], settings)

    test0_theta2 <- eval_util_L(settings, res0,
                                theta = array(runif(1 * I * J), c(1, I, J)),
                                cores = 1)
    checkmate::expect_data_frame(test0_theta2)
    expect_equal(colnames(test0_theta2), c(colnames(settings), "Utility"))
    expect_equal(test0_theta2[, -ncol(test0_theta2)], settings)

    test0_phi1 <- eval_util_L(settings, res0,
                              phi = matrix(runif(1 * I), 1, I),
                              cores = 1)
    checkmate::expect_data_frame(test0_phi1)
    expect_equal(colnames(test0_phi1), c(colnames(settings), "Utility"))
    expect_equal(test0_phi1[, -ncol(test0_phi1)], settings)

    test0_phi2 <- eval_util_L(settings, res0,
                              phi = array(runif(1 * I * J), c(1, I, J)),
                              cores = 1)
    checkmate::expect_data_frame(test0_phi2)
    expect_equal(colnames(test0_phi2), c(colnames(settings), "Utility"))
    expect_equal(test0_phi2[, -ncol(test0_phi2)], settings)
})

test_that("eval_util_R() outputs a data frame with the additional Utility column", {
    settings <- data.frame(J = rep(1, 3), K = rep(1, 3), N = rep(1, 3), x = NA)

    # Null model
    test0 <- eval_util_R(settings, res0, cores = 1)
    checkmate::expect_data_frame(test0)
    expect_equal(colnames(test0), c(colnames(settings), "Utility"))
    expect_equal(test0[, -ncol(test0)], settings)

    test0x <- eval_util_R(settings, res0, N_rep = 2, cores = 1)
    checkmate::expect_data_frame(test0x)
    expect_equal(colnames(test0x), c(colnames(settings), "Utility"))
    expect_equal(test0x[, -ncol(test0x)], settings)

    # Model with site covariates
    test1 <- eval_util_R(settings, res1, cores = 1)
    checkmate::expect_data_frame(test1)
    expect_equal(colnames(test1), c(colnames(settings), "Utility"))
    expect_equal(test1[, -ncol(test1)], settings)

    test1x <- eval_util_R(settings, res1, N_rep = 2, cores = 1)
    checkmate::expect_data_frame(test1x)
    expect_equal(colnames(test1x), c(colnames(settings), "Utility"))
    expect_equal(test1x[, -ncol(test1x)], settings)

    test1a <- eval_util_R(settings, res1a, cores = 1)
    checkmate::expect_data_frame(test1a)
    expect_equal(colnames(test1a), c(colnames(settings), "Utility"))
    expect_equal(test1a[, -ncol(test1a)], settings)

    test1ax <- eval_util_R(settings, res1a, N_rep = 2, cores = 1)
    checkmate::expect_data_frame(test1ax)
    expect_equal(colnames(test1ax), c(colnames(settings), "Utility"))
    expect_equal(test1ax[, -ncol(test1ax)], settings)

    test1b <- eval_util_R(settings, res1b, cores = 1)
    checkmate::expect_data_frame(test1b)
    expect_equal(colnames(test1b), c(colnames(settings), "Utility"))
    expect_equal(test1b[, -ncol(test1b)], settings)

    test1bx <- eval_util_R(settings, res1b, N_rep = 2, cores = 1)
    checkmate::expect_data_frame(test1bx)
    expect_equal(colnames(test1bx), c(colnames(settings), "Utility"))
    expect_equal(test1bx[, -ncol(test1bx)], settings)

    test2 <- eval_util_R(settings, res2, cores = 1)
    checkmate::expect_data_frame(test2)
    expect_equal(colnames(test2), c(colnames(settings), "Utility"))
    expect_equal(test2[, -ncol(test2)], settings)

    test2x <- eval_util_R(settings, res2, N_rep = 2, cores = 1)
    checkmate::expect_data_frame(test2x)
    expect_equal(colnames(test2x), c(colnames(settings), "Utility"))
    expect_equal(test2x[, -ncol(test2x)], settings)

    test3 <- eval_util_R(settings, res3, cores = 1)
    checkmate::expect_data_frame(test3)
    expect_equal(colnames(test3), c(colnames(settings), "Utility"))
    expect_equal(test3[, -ncol(test3)], settings)

    test3x <- eval_util_R(settings, res3, N_rep = 2, cores = 1)
    checkmate::expect_data_frame(test3x)
    expect_equal(colnames(test3x), c(colnames(settings), "Utility"))
    expect_equal(test3x[, -ncol(test3x)], settings)

    # Model with species covariates
    test4 <- eval_util_R(settings, res4, cores = 1)
    checkmate::expect_data_frame(test4)
    expect_equal(colnames(test4), c(colnames(settings), "Utility"))
    expect_equal(test4[, -ncol(test4)], settings)

    test4x <- eval_util_R(settings, res4, N_rep = 2, cores = 1)
    checkmate::expect_data_frame(test4x)
    expect_equal(colnames(test4x), c(colnames(settings), "Utility"))
    expect_equal(test4x[, -ncol(test4x)], settings)

    test5 <- eval_util_R(settings, res5, cores = 1)
    checkmate::expect_data_frame(test5)
    expect_equal(colnames(test5), c(colnames(settings), "Utility"))
    expect_equal(test5[, -ncol(test5)], settings)

    test5x <- eval_util_R(settings, res5, N_rep = 2, cores = 1)
    checkmate::expect_data_frame(test5x)
    expect_equal(colnames(test5x), c(colnames(settings), "Utility"))
    expect_equal(test5x[, -ncol(test5x)], settings)

    test6 <- eval_util_R(settings, res6, cores = 1)
    checkmate::expect_data_frame(test6)
    expect_equal(colnames(test6), c(colnames(settings), "Utility"))
    expect_equal(test6[, -ncol(test6)], settings)

    test6x <- eval_util_R(settings, res6, N_rep = 2, cores = 1)
    checkmate::expect_data_frame(test6x)
    expect_equal(colnames(test6x), c(colnames(settings), "Utility"))
    expect_equal(test6x[, -ncol(test6x)], settings)

    # Use psi, theta, phi arguments
    test0_psi1 <- eval_util_R(settings, res0,
                              psi = matrix(runif(1 * I), 1, I),
                              cores = 1)
    checkmate::expect_data_frame(test0_psi1)
    expect_equal(colnames(test0_psi1), c(colnames(settings), "Utility"))
    expect_equal(test0_psi1[, -ncol(test0_psi1)], settings)

    test0_psi2 <- eval_util_R(settings, res0,
                              psi = array(runif(1 * I * J), c(1, I, J)),
                              cores = 1)
    checkmate::expect_data_frame(test0_psi2)
    expect_equal(colnames(test0_psi2), c(colnames(settings), "Utility"))
    expect_equal(test0_psi2[, -ncol(test0_psi2)], settings)

    test0_theta1 <- eval_util_R(settings, res0,
                                theta = matrix(runif(1 * I), 1, I),
                                cores = 1)
    checkmate::expect_data_frame(test0_theta1)
    expect_equal(colnames(test0_theta1), c(colnames(settings), "Utility"))
    expect_equal(test0_theta1[, -ncol(test0_theta1)], settings)

    test0_theta2 <- eval_util_R(settings, res0,
                                theta = array(runif(1 * I * J), c(1, I, J)),
                                cores = 1)
    checkmate::expect_data_frame(test0_theta2)
    expect_equal(colnames(test0_theta2), c(colnames(settings), "Utility"))
    expect_equal(test0_theta2[, -ncol(test0_theta2)], settings)

    test0_phi1 <- eval_util_R(settings, res0,
                              phi = matrix(runif(1 * I), 1, I),
                              cores = 1)
    checkmate::expect_data_frame(test0_phi1)
    expect_equal(colnames(test0_phi1), c(colnames(settings), "Utility"))
    expect_equal(test0_phi1[, -ncol(test0_phi1)], settings)

    test0_phi2 <- eval_util_R(settings, res0,
                              phi = array(runif(1 * I * J), c(1, I, J)),
                              cores = 1)
    checkmate::expect_data_frame(test0_phi2)
    expect_equal(colnames(test0_phi2), c(colnames(settings), "Utility"))
    expect_equal(test0_phi2[, -ncol(test0_phi2)], settings)
})


### Tests for list_cond_L ------------------------------------------------------
I <- 2  # Number of species
J <- 50 # Number of sites
K <- 2  # Number of replicates
data <- occumbData(
    y = array(sample.int(I * J * K), dim = c(I, J, K)),
    spec_cov = list(cov1 = rnorm(I)),
    site_cov = list(cov2 = rnorm(J),
                    cov3 = factor(1:J)),
    repl_cov = list(cov4 = matrix(rnorm(J * K), J, K)))
res0 <- occumb(data = data,
               n.chains = 1, n.adapt = 0, n.burnin = 0,
               n.thin = 1, n.iter = 10, verbose = FALSE)
budget <- 850000; lambda1 <- 0.01; lambda2 <- 5000
test <- list_cond_L(budget, lambda1, lambda2, res0)

test_that("list_cond_L() outputs a data frame with correct columns", {
    checkmate::expect_data_frame(test)
    expect_equal(colnames(test), c("budget", "lambda1", "lambda2", "K", "N"))
})

test_that("Elements of list_cond_L() output are correct", {
    max_K <- floor(budget / (lambda2 * J))
    N <- (budget - lambda2 * J * seq_len(max_K)) / (lambda1 * J * seq_len(max_K))
    expect_equal(nrow(test), max_K)
    expect_equal(test$budget, rep(budget, max_K))
    expect_equal(test$lambda1, rep(lambda1, max_K))
    expect_equal(test$lambda2, rep(lambda2, max_K))
    expect_equal(test$K, seq_len(max_K))
    expect_equal(test$N, N)
})

test_that("K argument of list_cond_L() work correctly", {
    testK <- c(1, 3)
    N <- (budget - lambda2 * J * testK) / (lambda1 * J * testK)
    test <- list_cond_L(budget, lambda1, lambda2, res0, K = testK)
    expect_equal(nrow(test), length(testK))
    expect_equal(test$budget, rep(budget, length(testK)))
    expect_equal(test$lambda1, rep(lambda1, length(testK)))
    expect_equal(test$lambda2, rep(lambda2, length(testK)))
    expect_equal(test$K, testK)
    expect_equal(test$N, N)
})

test_that("Quality controls for list_cond_L() work correctly", {
    max_K <- floor(budget / (lambda2 * J))
    expect_error(list_cond_L(-1, lambda1, lambda2, res0),
                 "Negative 'budget' value.")
    expect_error(list_cond_L(budget, -1, lambda2, res0),
                 "Negative 'lambda1' value.")
    expect_error(list_cond_L(budget, lambda1, -1, res0),
                 "Negative 'lambda2' value.")
    expect_error(list_cond_L(budget, lambda1, lambda2, 0),
                 "An occumbFit class object is expected for 'fit'")
    expect_error(list_cond_L(0, lambda1, lambda2, res0),
                 "Impossible to have > 0 replicates per site under the given budget, cost, and the number of sites.")
    expect_error(list_cond_L(budget, lambda1, lambda2, res0, K = c(0, 1)),
                 "'K' contains values less than one.")
    expect_error(list_cond_L(budget, lambda1, lambda2, res0, K = seq_len(max_K + 1)),
                 paste("A value of 'K' greater than",
                       max_K,
                       "is not feasible under the given budget, cost, and the number of sites."))
})


### Tests for list_cond_R ------------------------------------------------------
budget <- 100000; lambda1 <- 0.01; lambda2 <- 5000; lambda3 <- 5000
max_K <- find_maxK(budget, lambda2, lambda3)
test <- list_cond_R(budget, lambda1, lambda2, lambda3)

test_that("list_cond_R() outputs a data frame with correct columns", {
    checkmate::expect_data_frame(test)
    expect_equal(colnames(test),
                 c("budget", "lambda1", "lambda2", "lambda3", "J", "K", "N"))
})

test_that("Elements of list_cond_R() output are correct", {
    J_valid <- list(); J <- seq_len(100)
    for (k in seq_len(max_K))
        J_valid[[k]] <- J[budget - lambda2 * J * k - lambda3 * J > 0]
    nrow_ans <- length(unlist(J_valid))
    J_ans <- K_ans <- vector()
    for (k in seq_len(max_K)) {
        J_ans <- c(J_ans, J_valid[[k]])
        K_ans <- c(K_ans, rep(k, sapply(J_valid, length)[k]))
    }
    N_ans <- (budget - lambda2 * J_ans * K_ans - lambda3 * J_ans) / (lambda1 * J_ans * K_ans)

    expect_equal(nrow(test), nrow_ans)
    expect_equal(test$budget, rep(budget, nrow_ans))
    expect_equal(test$lambda1, rep(lambda1, nrow_ans))
    expect_equal(test$lambda2, rep(lambda2, nrow_ans))
    expect_equal(test$lambda3, rep(lambda2, nrow_ans))
    expect_equal(test$J, J_ans)
    expect_equal(test$K, K_ans)
    expect_equal(test$N, N_ans)
})

test_that("J argument of list_cond_R() work correctly", {
    J <- seq(2, 10, 2)
    J_valid <- list()
    for (k in seq_len(max_K))
        J_valid[[k]] <- J[budget - lambda2 * J * k - lambda3 * J > 0]
    nrow_ans <- length(unlist(J_valid))
    J_ans <- K_ans <- vector()
    for (k in seq_len(max_K)) {
        J_ans <- c(J_ans, J_valid[[k]])
        K_ans <- c(K_ans, rep(k, sapply(J_valid, length)[k]))
    }
    N_ans <- (budget - lambda2 * J_ans * K_ans - lambda3 * J_ans) / (lambda1 * J_ans * K_ans)

    testJ <- list_cond_R(budget, lambda1, lambda2, lambda3, J = J)
    expect_equal(nrow(testJ), nrow_ans)
    expect_equal(testJ$budget, rep(budget, nrow_ans))
    expect_equal(testJ$lambda1, rep(lambda1, nrow_ans))
    expect_equal(testJ$lambda2, rep(lambda2, nrow_ans))
    expect_equal(testJ$lambda3, rep(lambda2, nrow_ans))
    expect_equal(testJ$J, J_ans)
    expect_equal(testJ$K, K_ans)
    expect_equal(testJ$N, N_ans)
})

test_that("K argument of list_cond_R() work correctly", {
    K <- c(1, 3)
    J_valid <- list(); J <- seq_len(100)
    for (k in K)
        J_valid[[k]] <- J[budget - lambda2 * J * k - lambda3 * J > 0]
    nrow_ans <- length(unlist(J_valid))
    J_ans <- K_ans <- vector()
    for (k in K) {
        J_ans <- c(J_ans, J_valid[[k]])
        K_ans <- c(K_ans, rep(k, sapply(J_valid, length)[k]))
    }
    N_ans <- (budget - lambda2 * J_ans * K_ans - lambda3 * J_ans) / (lambda1 * J_ans * K_ans)

    testK <- list_cond_R(budget, lambda1, lambda2, lambda3, K = K)
    expect_equal(nrow(testK), nrow_ans)
    expect_equal(testK$budget, rep(budget, nrow_ans))
    expect_equal(testK$lambda1, rep(lambda1, nrow_ans))
    expect_equal(testK$lambda2, rep(lambda2, nrow_ans))
    expect_equal(testK$lambda3, rep(lambda2, nrow_ans))
    expect_equal(testK$J, J_ans)
    expect_equal(testK$K, K_ans)
    expect_equal(testK$N, N_ans)
})

test_that("Quality controls for list_cond_R() work correctly", {
    expect_error(list_cond_R(-1, lambda1, lambda2, lambda3),
                 "Negative 'budget' value.")
    expect_error(list_cond_R(budget, -1, lambda2, lambda3),
                 "Negative 'lambda1' value.")
    expect_error(list_cond_R(budget, lambda1, -1, lambda3),
                 "Negative 'lambda2' value.")
    expect_error(list_cond_R(budget, lambda1, lambda2, -1),
                 "Negative 'lambda3' value.")
    expect_error(list_cond_R(budget, lambda1, lambda2, lambda3, J = c(0, 1)),
                 "'J' contains values less than one.")
    expect_error(list_cond_R(budget, lambda1, lambda2, lambda3, K = c(0, 1)),
                 "'K' contains values less than one.")
    expect_error(list_cond_R(budget, lambda1, lambda2, lambda3, K = c(2, 1)),
                 "'K' must be in ascending order.")
    expect_error(list_cond_R(budget, lambda1, lambda2, lambda3, K = c(20)),
                 paste("No valid combination of 'J' and 'K' under the given budget and cost."))
})


### Tests for find_max_J/K -----------------------------------------------------
test_that("find_maxJ/K() returns correct value", {
    budget <- 100000; lambda2 <- 5000; lambda3 <- 5000
    J <- seq_len(1E3)
    J_ans <- max(J[budget - lambda2 * J - lambda3 * J > 0])
    K <- seq_len(1E3)
    K_ans <- max(K[budget - lambda2 * K - lambda3 > 0])

    expect_equal(find_maxJ(budget, lambda2, lambda3), J_ans)
    expect_equal(find_maxK(budget, lambda2, lambda3), K_ans)
})

test_that("find_maxJ/K() returns zero when budget is too small", {
    expect_equal(find_maxJ(10, lambda2, lambda3), 0)
    expect_equal(find_maxK(10, lambda2, lambda3), 0)
})

test_that("find_maxJ/K() throws an error when the budget is too large", {
    expect_error(find_maxJ(1E16, lambda2, lambda3),
                 "Maximum `J` value seems too large under the specified budget and cost values: consider using the `J` argument to specify a smaller set of `J` values of interest.")
    expect_error(find_maxK(1E16, lambda2, lambda3),
                 "Maximum `K` value seems too large under the specified budget and cost values: consider using the `K` argument to specify a smaller set of `K` values of interest.")
})


### Tests for check_args_eval_util_L -------------------------------------------
I <- dim(res0@data@y)[1]
J <- dim(res0@data@y)[2]
df_test  <- data.frame(K = rep(1, 2), N = rep(1, 2))
arr_test <- array(1, dim = c(1, I, J))
mat_test <- array(1, dim = c(1, I))

test_that("check_args_eval_util_L() blocks inappropriate settings", {
    expect_error(check_args_eval_util_L(
                    data.frame(Kx = rep(1, 2), N = rep(1, 2)), res0,
                    NULL, NULL, NULL),
                 "The 'settings' argument does not contain column 'K'.")
    expect_error(check_args_eval_util_L(
                    data.frame(K = rep(1, 2), Nx = rep(1, 2)), res0,
                    NULL, NULL, NULL),
                 "The 'settings' argument does not contain column 'N'.")
    expect_error(check_args_eval_util_L(
                    data.frame(K = rep(0, 2), N = rep(1, 2)), res0,
                    NULL, NULL, NULL),
                 "'K' contains values less than one.")
    expect_error(check_args_eval_util_L(
                    data.frame(K = rep(1, 2), N = rep(0, 2)), res0,
                    NULL, NULL, NULL),
                 "'N' contains values less than one.")
})

test_that("check_args_eval_util_L() allows sufficient arguments", {
    expect_invisible(
        check_args_eval_util_L(df_test, res0, NULL, NULL, NULL))
    expect_invisible(
        check_args_eval_util_L(df_test, NULL, arr_test, arr_test, arr_test))
    expect_invisible(
        check_args_eval_util_L(df_test, NULL, arr_test, mat_test, arr_test))
    expect_invisible(
        check_args_eval_util_L(df_test, NULL, arr_test, arr_test, mat_test))
    expect_invisible(
        check_args_eval_util_L(df_test, NULL, arr_test, mat_test, mat_test))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, arr_test, NULL, NULL))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, NULL, arr_test, NULL))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, NULL, mat_test, NULL))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, NULL, NULL, arr_test))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, NULL, NULL, mat_test))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, arr_test, arr_test, NULL))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, arr_test, mat_test, NULL))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, arr_test, NULL, arr_test))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, arr_test, NULL, mat_test))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, NULL, arr_test, arr_test))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, NULL, mat_test, arr_test))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, NULL, arr_test, mat_test))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, NULL, mat_test, mat_test))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, arr_test, arr_test, arr_test))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, arr_test, mat_test, arr_test))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, arr_test, arr_test, mat_test))
    expect_invisible(
        check_args_eval_util_L(df_test, res0, arr_test, mat_test, mat_test))
})

test_that("check_args_eval_util_L() blocks insufficient arguments", {
    error_message <- "Parameter values are not fully specified: use fit argument or otherwise use all of z, theta, phi arguments."
    expect_error(check_args_eval_util_L(
                    df_test, NULL, NULL, NULL, NULL),
                 error_message)
    expect_error(check_args_eval_util_L(
                    df_test, NULL, arr_test, NULL, NULL),
                 error_message)
    expect_error(check_args_eval_util_L(
                    df_test, NULL, NULL, arr_test, NULL),
                 error_message)
    expect_error(check_args_eval_util_L(
                    df_test, NULL, NULL, mat_test, NULL),
                 error_message)
    expect_error(check_args_eval_util_L(
                    df_test, NULL, NULL, NULL, arr_test),
                 error_message)
    expect_error(check_args_eval_util_L(
                    df_test, NULL, NULL, NULL, mat_test),
                 error_message)
    expect_error(check_args_eval_util_L(
                    df_test, NULL, NULL, arr_test, arr_test),
                 error_message)
    expect_error(check_args_eval_util_L(
                    df_test, NULL, NULL, mat_test, arr_test),
                 error_message)
    expect_error(check_args_eval_util_L(
                    df_test, NULL, NULL, arr_test, mat_test),
                 error_message)
    expect_error(check_args_eval_util_L(
                    df_test, NULL, NULL, mat_test, mat_test),
                 error_message)
    expect_error(check_args_eval_util_L(
                    df_test, NULL, arr_test, NULL, arr_test),
                 error_message)
    expect_error(check_args_eval_util_L(
                    df_test, NULL, arr_test, NULL, mat_test),
                 error_message)
    expect_error(check_args_eval_util_L(
                    df_test, NULL, arr_test, arr_test, NULL),
                 error_message)
    expect_error(check_args_eval_util_L(
                    df_test, NULL, arr_test, mat_test, NULL),
                 error_message)
})

test_that("check_args_eval_util_L() blocks models with replicate-specific parameters", {
    # Replicate-specific model is supplied, but not phi/theta
    expect_error(check_args_eval_util_L(df_test, res7, NULL, NULL, NULL),
                 "'fit' contains replicate-specific phi: specify appropriate phi values via the 'phi' argument to run.")
    expect_error(check_args_eval_util_L(df_test, res8, NULL, NULL, NULL),
                 "'fit' contains replicate-specific theta: specify appropriate theta values via the 'theta' argument to run.")

    # Both replicate-specific phi model and phi samples are supplied
    I <- dim(res7@data@y)[1]
    J <- dim(res7@data@y)[2]
    mat_test <- array(1, dim = c(1, I))
    arr_test <- array(1, dim = c(1, I, J))
    expect_no_error(check_args_eval_util_L(df_test, res7, NULL, NULL, mat_test))
    expect_no_error(check_args_eval_util_L(df_test, res7, NULL, NULL, arr_test))

    # Both replicate-specific theta model and theta samples are supplied
    I <- dim(res8@data@y)[1]
    J <- dim(res8@data@y)[2]
    mat_test <- array(1, dim = c(1, I))
    arr_test <- array(1, dim = c(1, I, J))
    expect_no_error(check_args_eval_util_L(df_test, res8, NULL, mat_test, NULL))
    expect_no_error(check_args_eval_util_L(df_test, res8, NULL, arr_test, NULL))
})

test_that("check_args_eval_util_L() blocks z with sites occupied by no species", {
    expect_error(check_args_eval_util_L(
                    df_test, res0, array(0, dim = c(1, I, J)), NULL, NULL),
                 "The given 'z' array contains case\\(s\\) where no species occupy a site; see, for example, that z\\[1, , 1\\] is a zero vector")
})

test_that("check_args_eval_util_L() blocks dimension mismatch between z, theta, phi and fit", {
    expect_error(check_args_eval_util_L(
                    df_test, res0, array(1, dim = c(1, I + 1, J)), NULL, NULL),
                 paste0("Mismatch in species dimension: dim\\(z\\)\\[2\\] must be ", I, ".\n"))
    expect_error(check_args_eval_util_L(
                    df_test, res0, array(1, dim = c(1, I, J + 1)), NULL, NULL),
                 paste0("Mismatch in site dimension: dim\\(z\\)\\[3\\] must be ", J, ".\n"))

    expect_error(check_args_eval_util_L(
                    df_test, res0, NULL, array(1, dim = c(1, I + 1, J)), NULL),
                 paste0("Mismatch in species dimension: dim\\(theta\\)\\[2\\] must be ", I, ".\n"))
    expect_error(check_args_eval_util_L(
                    df_test, res0, NULL, array(1, dim = c(1, I, J + 1)), NULL),
                 paste0("Mismatch in site dimension: dim\\(theta\\)\\[3\\] must be ", J, ".\n"))

    expect_error(check_args_eval_util_L(
                    df_test, res0, NULL, NULL, array(1, dim = c(1, I + 1, J))),
                 paste0("Mismatch in species dimension: dim\\(phi\\)\\[2\\] must be ", I, ".\n"))
    expect_error(check_args_eval_util_L(
                    df_test, res0, NULL, NULL, array(1, dim = c(1, I, J + 1))),
                 paste0("Mismatch in site dimension: dim\\(phi\\)\\[3\\] must be ", J, ".\n"))

    expect_error(check_args_eval_util_L(df_test, NULL,
                    arr_test,
                    array(1, dim = c(1, I + 1, J)),
                    arr_test),
                 "Mismatch in species dimension: dim\\(z\\)\\[2\\], dim\\(theta\\)\\[2\\], and dim\\(phi\\)\\[2\\] must be equal.")
    expect_error(check_args_eval_util_L(df_test, NULL,
                    arr_test,
                    array(1, dim = c(1, I + 1, J)),
                    mat_test),
                 "Mismatch in species dimension: dim\\(z\\)\\[2\\], dim\\(theta\\)\\[2\\], and dim\\(phi\\)\\[2\\] must be equal.")
    expect_error(check_args_eval_util_L(df_test, NULL,
                    arr_test,
                    arr_test,
                    array(1, dim = c(1, I + 1, J))),
                 "Mismatch in species dimension: dim\\(z\\)\\[2\\], dim\\(theta\\)\\[2\\], and dim\\(phi\\)\\[2\\] must be equal.")
    expect_error(check_args_eval_util_L(df_test, NULL,
                    arr_test,
                    mat_test,
                    array(1, dim = c(1, I + 1, J))),
                 "Mismatch in species dimension: dim\\(z\\)\\[2\\], dim\\(theta\\)\\[2\\], and dim\\(phi\\)\\[2\\] must be equal.")

    expect_error(check_args_eval_util_L(df_test, NULL,
                    arr_test,
                    array(1, dim = c(1, I, J + 1)),
                    arr_test),
                 "Mismatch in site dimension: dim\\(z\\)\\[3\\], dim\\(theta\\)\\[3\\], and dim\\(phi\\)\\[3\\] must be equal.")
    expect_error(check_args_eval_util_L(df_test, NULL,
                    arr_test,
                    array(1, dim = c(1, I, J + 1)),
                    mat_test),
                 "Mismatch in site dimension: dim\\(z\\)\\[3\\] and dim\\(theta\\)\\[3\\] must be equal.")
    expect_error(check_args_eval_util_L(df_test, NULL,
                    arr_test,
                    arr_test,
                    array(1, dim = c(1, I, J + 1))),
                 "Mismatch in site dimension: dim\\(z\\)\\[3\\], dim\\(theta\\)\\[3\\], and dim\\(phi\\)\\[3\\] must be equal.")
    expect_error(check_args_eval_util_L(df_test, NULL,
                    arr_test,
                    mat_test,
                    array(1, dim = c(1, I, J + 1))),
                 "Mismatch in site dimension: dim\\(z\\)\\[3\\] and dim\\(phi\\)\\[3\\] must be equal.")

    expect_error(check_args_eval_util_L(df_test, NULL,
                    arr_test,
                    array(1, dim = c(1, I, J + 1)),
                    array(1, dim = c(1, I))),
                 "Mismatch in site dimension: dim\\(z\\)\\[3\\] and dim\\(theta\\)\\[3\\] must be equal.")
    expect_error(check_args_eval_util_L(df_test, NULL,
                    arr_test,
                    array(1, dim = c(1, I)),
                    array(1, dim = c(1, I, J + 1))),
                 "Mismatch in site dimension: dim\\(z\\)\\[3\\] and dim\\(phi\\)\\[3\\] must be equal.")
})


### Tests for check_args_eval_util_R ------------------------------------------
I <- dim(res0@data@y)[1]
J <- dim(res0@data@y)[2]
df_test  <- data.frame(J = rep(1, 2), K = rep(1, 2), N = rep(1, 2))
arr_test <- array(1, dim = c(1, I, J))
mat_test <- array(1, dim = c(1, I))

test_that("check_args_eval_util_R() blocks inappropriate settings", {
    expect_error(check_args_eval_util_R(
                    data.frame(Jx = rep(1, 2), K = rep(1, 2), N = rep(1, 2)),
                    res0),
                 "The 'settings' argument does not contain column 'J'.")
    expect_error(check_args_eval_util_R(
                    data.frame(J = rep(1, 2), Kx = rep(1, 2), N = rep(1, 2)),
                    res0),
                 "The 'settings' argument does not contain column 'K'.")
    expect_error(check_args_eval_util_R(
                    data.frame(J = rep(1, 2),
                               K = rep(1, 2),
                               Nx = rep(1, 2)),
                    res0),
                 "The 'settings' argument does not contain column 'N'.")
    expect_error(check_args_eval_util_R(
                    data.frame(J = rep(0, 2),
                               K = rep(1, 2),
                               N = rep(1, 2)),
                    res0),
                 "'J' contains values less than one.")
    expect_error(check_args_eval_util_R(
                    data.frame(J = rep(1, 2),
                               K = rep(0, 2),
                               N = rep(1, 2)),
                    res0),
                 "'K' contains values less than one.")
    expect_error(check_args_eval_util_R(
                    data.frame(J = rep(1, 2),
                               K = rep(1, 2),
                               N = rep(0, 2)),
                    res0),
                 "'N' contains values less than one.")
})

test_that("check_args_eval_util_R() allows sufficient arguments", {
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0, NULL, NULL, NULL))
    expect_invisible(check_args_eval_util_R(df_test,
                                            NULL,
                                            arr_test,
                                            arr_test,
                                            arr_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            NULL,
                                            mat_test,
                                            arr_test,
                                            arr_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            NULL,
                                            arr_test,
                                            mat_test,
                                            arr_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            NULL,
                                            arr_test,
                                            arr_test,
                                            mat_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            NULL,
                                            mat_test,
                                            mat_test,
                                            arr_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            NULL,
                                            mat_test,
                                            arr_test,
                                            mat_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            NULL,
                                            arr_test,
                                            mat_test,
                                            mat_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            NULL,
                                            mat_test,
                                            mat_test,
                                            mat_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            arr_test,
                                            NULL,
                                            NULL))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            mat_test,
                                            NULL,
                                            NULL))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            NULL,
                                            arr_test,
                                            NULL))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            NULL,
                                            mat_test,
                                            NULL))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            NULL,
                                            NULL,
                                            arr_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            NULL,
                                            NULL,
                                            mat_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            arr_test,
                                            arr_test,
                                            NULL))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            mat_test,
                                            arr_test,
                                            NULL))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            arr_test,
                                            mat_test,
                                            NULL))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            mat_test,
                                            mat_test,
                                            NULL))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            arr_test,
                                            NULL,
                                            arr_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            mat_test,
                                            NULL,
                                            arr_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            arr_test,
                                            NULL,
                                            mat_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            mat_test,
                                            NULL,
                                            mat_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            NULL,
                                            arr_test,
                                            arr_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            NULL,
                                            mat_test,
                                            arr_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            NULL,
                                            arr_test,
                                            mat_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            NULL,
                                            mat_test,
                                            mat_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            arr_test,
                                            arr_test,
                                            arr_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            mat_test,
                                            arr_test,
                                            arr_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            arr_test,
                                            mat_test,
                                            arr_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            arr_test,
                                            arr_test,
                                            mat_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            mat_test,
                                            mat_test,
                                            arr_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            mat_test,
                                            arr_test,
                                            mat_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            arr_test,
                                            mat_test,
                                            mat_test))
    expect_invisible(check_args_eval_util_R(df_test,
                                            res0,
                                            mat_test,
                                            mat_test,
                                            mat_test))
})

test_that("check_args_eval_util_R() blocks insufficient arguments", {
    error_message <- "Parameter values are not fully specified: use fit argument or otherwise use all of psi, theta, phi arguments."
    expect_error(check_args_eval_util_R(df_test, NULL, NULL, NULL, NULL),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        NULL,
                                        NULL),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        NULL,
                                        NULL),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        NULL,
                                        arr_test,
                                        NULL),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        NULL,
                                        mat_test,
                                        NULL),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        NULL,
                                        NULL,
                                        arr_test),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        NULL,
                                        NULL,
                                        mat_test),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        NULL,
                                        arr_test,
                                        arr_test),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        NULL,
                                        mat_test,
                                        arr_test),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        NULL,
                                        arr_test,
                                        mat_test),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        NULL,
                                        mat_test,
                                        mat_test),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        NULL,
                                        arr_test),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        NULL,
                                        arr_test),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        NULL,
                                        mat_test),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        NULL,
                                        mat_test),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        arr_test,
                                        NULL),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        arr_test,
                                        NULL),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        mat_test,
                                        NULL),
                 error_message)
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        mat_test,
                                        NULL),
                 error_message)
})

test_that("check_args_eval_util_R() blocks models with replicate-specific parameters", {
    # Replicate-specific model is supplied, but not phi/theta
    expect_error(check_args_eval_util_R(df_test, res7, NULL, NULL, NULL),
                 "'fit' contains replicate-specific phi: specify appropriate phi values via the 'phi' argument to run.")
    expect_error(check_args_eval_util_R(df_test, res8, NULL, NULL, NULL),
                 "'fit' contains replicate-specific theta: specify appropriate theta values via the 'theta' argument to run.")

    # Both replicate-specific phi model and phi samples are supplied
    I <- dim(res7@data@y)[1]
    J <- dim(res7@data@y)[2]
    mat_test <- array(1, dim = c(1, I))
    arr_test <- array(1, dim = c(1, I, J))
    expect_no_error(check_args_eval_util_R(df_test, res7, NULL, NULL, mat_test))
    expect_no_error(check_args_eval_util_R(df_test, res7, NULL, NULL, arr_test))

    # Both replicate-specific theta model and theta samples are supplied
    I <- dim(res8@data@y)[1]
    J <- dim(res8@data@y)[2]
    mat_test <- array(1, dim = c(1, I))
    arr_test <- array(1, dim = c(1, I, J))
    expect_no_error(check_args_eval_util_R(df_test, res8, NULL, mat_test, NULL))
    expect_no_error(check_args_eval_util_R(df_test, res8, NULL, arr_test, NULL))
})

test_that("check_args_eval_util_R() blocks species/site dimension mismatch between psi, theta, phi, and fit", {
    expect_error(check_args_eval_util_R(
                    df_test, res0, array(1, dim = c(1, I + 1, J)), NULL, NULL),
                 paste0("Mismatch in species dimension: dim\\(psi\\)\\[2\\] must be ", I, ".\n"))
    expect_error(check_args_eval_util_R(
                    df_test, res0, array(1, dim = c(1, I, J + 1)), NULL, NULL),
                 paste0("Mismatch in site dimension: dim\\(psi\\)\\[3\\] must be ", J, ".\n"))

    expect_error(check_args_eval_util_R(
                    df_test, res0, NULL, array(1, dim = c(1, I + 1, J)), NULL),
                 paste0("Mismatch in species dimension: dim\\(theta\\)\\[2\\] must be ", I, ".\n"))
    expect_error(check_args_eval_util_R(
                    df_test, res0, NULL, array(1, dim = c(1, I, J + 1)), NULL),
                 paste0("Mismatch in site dimension: dim\\(theta\\)\\[3\\] must be ", J, ".\n"))

    expect_error(check_args_eval_util_R(
                    df_test, res0, NULL, NULL, array(1, dim = c(1, I + 1, J))),
                 paste0("Mismatch in species dimension: dim\\(phi\\)\\[2\\] must be ", I, ".\n"))
    expect_error(check_args_eval_util_R(
                    df_test, res0, NULL, NULL, array(1, dim = c(1, I, J + 1))),
                 paste0("Mismatch in site dimension: dim\\(phi\\)\\[3\\] must be ", J, ".\n"))

    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        array(1, dim = c(1, I + 1, J)),
                                        arr_test),
                 "Mismatch in species dimension: dim\\(psi\\)\\[2\\], dim\\(theta\\)\\[2\\], and dim\\(phi\\)\\[2\\] must be equal.")
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        array(1, dim = c(1, I + 1, J)),
                                        arr_test),
                 "Mismatch in species dimension: dim\\(psi\\)\\[2\\], dim\\(theta\\)\\[2\\], and dim\\(phi\\)\\[2\\] must be equal.")
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        array(1, dim = c(1, I + 1, J)),
                                        mat_test),
                 "Mismatch in species dimension: dim\\(psi\\)\\[2\\], dim\\(theta\\)\\[2\\], and dim\\(phi\\)\\[2\\] must be equal.")
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        array(1, dim = c(1, I + 1, J)),
                                        mat_test),
                 "Mismatch in species dimension: dim\\(psi\\)\\[2\\], dim\\(theta\\)\\[2\\], and dim\\(phi\\)\\[2\\] must be equal.")
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        arr_test,
                                        array(1, dim = c(1, I + 1, J))),
                 "Mismatch in species dimension: dim\\(psi\\)\\[2\\], dim\\(theta\\)\\[2\\], and dim\\(phi\\)\\[2\\] must be equal.")
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        arr_test,
                                        array(1, dim = c(1, I + 1, J))),
                 "Mismatch in species dimension: dim\\(psi\\)\\[2\\], dim\\(theta\\)\\[2\\], and dim\\(phi\\)\\[2\\] must be equal.")
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        mat_test,
                                        array(1, dim = c(1, I + 1, J))),
                 "Mismatch in species dimension: dim\\(psi\\)\\[2\\], dim\\(theta\\)\\[2\\], and dim\\(phi\\)\\[2\\] must be equal.")
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        mat_test,
                                        array(1, dim = c(1, I + 1, J))),
                 "Mismatch in species dimension: dim\\(psi\\)\\[2\\], dim\\(theta\\)\\[2\\], and dim\\(phi\\)\\[2\\] must be equal.")

    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        array(1, dim = c(1, I, J + 1)),
                                        arr_test),
                 "Mismatch in site dimension: dim\\(psi\\)\\[3\\], dim\\(theta\\)\\[3\\], and dim\\(phi\\)\\[3\\] must be equal.")
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        array(1, dim = c(1, I, J + 1)),
                                        arr_test),
                 "Mismatch in site dimension: dim\\(theta\\)\\[3\\] and dim\\(phi\\)\\[3\\] must be equal.")
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        array(1, dim = c(1, I, J + 1)),
                                        mat_test),
                 "Mismatch in site dimension: dim\\(psi\\)\\[3\\] and dim\\(theta\\)\\[3\\] must be equal.")
    expect_no_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        array(1, dim = c(1, I, J + 1)),
                                        mat_test))
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        arr_test,
                                        array(1, dim = c(1, I, J + 1))),
                 "Mismatch in site dimension: dim\\(psi\\)\\[3\\], dim\\(theta\\)\\[3\\], and dim\\(phi\\)\\[3\\] must be equal.")
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        arr_test,
                                        array(1, dim = c(1, I, J + 1))),
                 "Mismatch in site dimension: dim\\(theta\\)\\[3\\] and dim\\(phi\\)\\[3\\] must be equal.")
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        mat_test,
                                        array(1, dim = c(1, I, J + 1))),
                 "Mismatch in site dimension: dim\\(psi\\)\\[3\\] and dim\\(phi\\)\\[3\\] must be equal.")
    expect_no_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        mat_test,
                                        array(1, dim = c(1, I, J + 1))))

    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        array(1, dim = c(1, I, J + 1)),
                                        array(1, dim = c(1, I))),
                 "Mismatch in site dimension: dim\\(psi\\)\\[3\\] and dim\\(theta\\)\\[3\\] must be equal.")
    expect_no_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        array(1, dim = c(1, I, J + 1)),
                                        array(1, dim = c(1, I))))
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        arr_test,
                                        array(1, dim = c(1, I)),
                                        array(1, dim = c(1, I, J + 1))),
                 "Mismatch in site dimension: dim\\(psi\\)\\[3\\] and dim\\(phi\\)\\[3\\] must be equal.")
    expect_no_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        mat_test,
                                        array(1, dim = c(1, I)),
                                        array(1, dim = c(1, I, J + 1))))
    expect_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        array(1, dim = c(1, I)),
                                        arr_test,
                                        array(1, dim = c(1, I, J + 1))),
                 "Mismatch in site dimension: dim\\(theta\\)\\[3\\] and dim\\(phi\\)\\[3\\] must be equal.")
    expect_no_error(check_args_eval_util_R(df_test,
                                        NULL,
                                        array(1, dim = c(1, I)),
                                        mat_test,
                                        array(1, dim = c(1, I, J + 1))))
})




### Tests for eutil ------------------------------------------------------------

I <- 20; J <- 5; K <- 4; N <- 100; M <- 2
seed  <- rnorm(1)
z     <- array(rbinom(M * I * J, 1, 0.8), dim = c(M, I, J))
theta <- array(runif(M * I * J, min = 0.8), dim = c(M, I, J))
phi   <- array(rgamma(M * I * J, 1), dim = c(M, I, J))

test_that("eutil() works as expected for local scale", {
    ## * Tests are available only for non-parallel computations *
    # N_rep = 1
    ans <- with_seed(seed, {
        util_rep <- vector(length = M)
        for (n in seq_len(M))
            util_rep[n] <-
                cutil_local(z[n, , ], theta[n, , ], phi[n, , ], K, N)
        mean(util_rep)}
    )
    expect_equal(
        with_seed(seed, eutil(z, theta, phi, K, N, scale = "local",
                              N_rep = 1, cores = 1)),
        ans
    )
    # rep > 1
    ans <- with_seed(seed, {
        util_rep <- vector(length = M * 2)
        n <- rep(seq_len(M), each = 2)
        for (m in seq_along(n))
            util_rep[m] <-
                cutil_local(z[n[m], , ], theta[n[m], , ], phi[n[m], , ], K, N)
        mean(util_rep)}
    )
    expect_equal(
        with_seed(seed, eutil(z, theta, phi, K, N, scale = "local",
                              N_rep = 2, cores = 1)),
        ans
    )
})

test_that("eutil() works as expected for regional scale", {
    ## * Tests are available only for non-parallel computations *
    # rep = 1
    ans <- with_seed(seed, {
        util_rep <- vector(length = M)
        for (n in seq_len(M))
            util_rep[n] <-
                cutil_regional(z[n, , ], theta[n, , ], phi[n, , ], K, N)
        mean(util_rep)}
    )
    expect_equal(
        with_seed(seed, eutil(z, theta, phi, K, N, scale = "regional",
                              N_rep = 1, cores = 1)),
        ans
    )
    # rep > 1
    ans <- with_seed(seed, {
        util_rep <- vector(length = M * 2)
        n <- rep(seq_len(M), each = 2)
        for (m in seq_along(n))
            util_rep[m] <-
                cutil_regional(z[n[m], , ], theta[n[m], , ], phi[n[m], , ], K, N)
        mean(util_rep)}
    )
    expect_equal(
        with_seed(seed, eutil(z, theta, phi, K, N, scale = "regional",
                              N_rep = 2, cores = 1)),
        ans
    )
})

test_that("eutil() runs in parallel", {
    expect_visible(
        eutil(z, theta, phi, K, N, scale = "local", N_rep = 1, cores = 2)
    )
    expect_visible(
        eutil(z, theta, phi, K, N, scale = "regional", N_rep = 1, cores = 2)
    )
})

test_that("eutil() outputs errors when executed in parallel", {
    theta <- array(1E-10, dim = c(M, I, J))
    suppressWarnings({
    expect_error(
        eutil(z, theta, phi, K, N, scale = "local", N_rep = 1, cores = 2)
    )
    })
    suppressWarnings({
    expect_error(
        eutil(z, theta, phi, K, N, scale = "regional", N_rep = 1, cores = 2)
    )
    })
})

test_that("eutil() warns when utility when some util_rep values could not be calculated", {
    phi   <- array(1E-10, dim = c(M, I, J))
    expect_warning(
        eutil(z, theta, phi, K, N, scale = "local", N_rep = 1, cores = 1),
        "Case\\(s\\) arose in the replicated simulation where 'Utility' could not be calculated and were ignored. This result may sometimes occur stochastically; try repeat running to see if the same warning occurs. If the same result occurs frequently, the given 'theta' or 'phi' values might need to be higher."
    )
    expect_warning(
        eutil(z, theta, phi, K, N, scale = "regional", N_rep = 1, cores = 1),
        "Case\\(s\\) arose in the replicated simulation where 'Utility' could not be calculated and were ignored. This result may sometimes occur stochastically; try repeat running to see if the same warning occurs. If the same result occurs frequently, the given 'theta' or 'phi' values might need to be higher."
    )
})

### Tests for cutil ------------------------------------------------------------

I <- 20; J <- 5; K <- 4; N <- 100
u <- r <- pi <- array(dim = c(I, J, K)); pi[1] <- NaN
while(any(is.nan(pi))) {
    z     <- matrix(rbinom(I * J, 1, 0.8), I, J)
    theta <- matrix(runif(I * J, min = 0.8), I, J)
    phi   <- matrix(rgamma(I * J, 1), I, J)
    seed  <- rnorm(1)
    set.seed(seed)

    for (k in seq_len(K)) {
        u[, , k] <- rbinom(I * J, 1, z * theta)
        r[, , k] <- rgamma(I * J, phi, 1)
    }
    for (k in seq_len(K)) {
        for (j in seq_len(J)) {
            for (i in seq_len(I)) {
                pi[i, j, k] <- u[i, j, k] * r[i, j, k] / sum(u[, j, k] * r[, j, k])
            }
        }
    }
}

test_that("cutil_local() works as expected", {
    ans <- with_seed(seed,
        sum(predict_detect_probs_local(predict_pi(z, theta, phi, K), N)) / J)
    expect_equal(with_seed(seed, cutil_local(z, theta, phi, K, N)), ans)
})

test_that("cutil_regional() works as expected", {
    ans <- with_seed(seed,
        sum(predict_detect_probs_regional(predict_pi(z, theta, phi, K), N)))
    expect_equal(with_seed(seed, cutil_regional(z, theta, phi, K, N)), ans)
})

test_that("predict_detect_probs_local() works as expected", {
    ans <- matrix(nrow = I, ncol = J)
    for (i in seq_len(I)) {
        for (j in seq_len(J)) {
            ans[i, j] <- 1 - prod((1 - pi[i, j, ])^N)
        }
    }
    expect_equal(predict_detect_probs_local(pi, N), ans)
})

test_that("predict_detect_probs_regional() works as expected", {
    ans <- vector(length = I)
    for (i in seq_len(I)) ans[i] <- 1 - prod((1 - pi[i, , ])^N)
    expect_equal(predict_detect_probs_regional(pi, N), ans)
})

### Tests for sample_z ---------------------------------------------------------
test_that("sample_z() stops when z sampling fails repeatedly", {
    expect_error(sample_z(rep(1E-10, 2)),
                 "Failed to generate valid 'z' values under the given parameter set. Providing 'psi' containing higher psi values may fix the issue.")
})

### Tests for sample_u ---------------------------------------------------------
test_that("sample_u() stops when u sampling fails repeatedly", {
    expect_error(sample_u(matrix(rep(1E-10, 2), 2, 1)),
                 "Failed to generate valid 'u' values under the given parameter set. Providing 'theta' containing higher theta values may fix the issue.")
})

