### Test for validate_occumbData() ---------------------------------------------
test_that("Dimension check for y works", {
    expect_error(new("occumbData", y = array(1:4, dim = rep(2, 1))),
                 "'y' should be a 3D-array.")
    expect_error(new("occumbData", y = array(1:4, dim = rep(2, 2))),
                 "'y' should be a 3D-array.")
    expect_error(new("occumbData", y = array(1:16, dim = rep(2, 4))),
                 "'y' should be a 3D-array.")
})

test_that("Check for missing values for y works", {
    expect_error(new("occumbData", y = array(c(NA, 1:7), dim = rep(2, 3))),
                 "'y' contains missing value.")
})

test_that("Integer check works", {
    expect_error(new("occumbData", y = array(1:8 + 0.1, dim = rep(2, 3))),
                 "'y' contains non-integer value.")
})

test_that("Check for covariate name overlap works", {
    expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                                   spec_cov = list(a = NULL, b = NULL),
                                   site_cov = list(a = NULL)),
                 "Duplicated covariate names are not allowed: 'a'")
})

test_that("Dimension check for covariates works", {
    expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                                   spec_cov = list(a = rep(1, 1))),
                 "Length of 'a' should match the number of species.")
    expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                                   site_cov = list(b = rep(1, 3))),
                 "Length of 'b' should match the number of sites.")
    expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                                   repl_cov = list(c = rep(1, 3))),
                 "'c' should be a matrix with J rows and K columns.")
    expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                                   repl_cov = list(c = matrix(1:2, 1, 2))),
                 "'c' should be a matrix with J rows and K columns.")
})

test_that("Check for covariate missing values works", {
    expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                                   spec_cov = list(a = c(1, NA))),
                 "'spec_cov' contains missing value.")
    expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                                   site_cov = list(b = c(1, NA))),
                 "'site_cov' contains missing value.")
    expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                                   repl_cov = list(c = matrix(c(1:3, NA), 2, 2))),
                 "'repl_cov' contains missing value.")
})

### Test for check_covariate_mode() --------------------------------------------
test_that("Check for unacceptable covariate mode works", {
    spec_cov <- list(a = rnorm(1), b = c(0 + 1i))
    site_cov <- list(c = rnorm(1), d = c(0 + 1i))
    repl_cov <- list(e = rnorm(1), f = c(0 + 1i))
    expect_error(check_covariate_mode(spec_cov, site_cov, repl_cov),
                 sprintf("Unacceptable mode: the following covariates must be numeric, factor, or character. \n %s",
                         paste(sprintf("%s: %s", c("b", "d", "f"), rep("complex", 3)),
                               collapse = "; ")))
})

