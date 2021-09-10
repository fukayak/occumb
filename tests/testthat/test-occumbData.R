### Test for validate_occumbData() --------------------------------------------
test_that("Dimension check for y works", {
    expect_error(new("occumbData", y = array(1:4, dim = rep(2, 2))),
                 "'y' should be a 3D-array.")
    expect_error(new("occumbData", y = array(1:16, dim = rep(2, 4))),
                 "'y' should be a 3D-array.")
})

test_that("Check for missing values for y works", {
    expect_error(new("occumbData", y = array(c(NA, 1:7), dim = rep(2, 3))),
                 "Missing values are not allowed in 'y'.")
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
                 "Length of 'c' should match the number of replicates.")
})

test_that("Check for covariate missing values works", {
    expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                                   spec_cov = list(a = c(1, NA))),
                 "Missing values are not allowed in 'spec_cov'.")
    expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                                   site_cov = list(b = c(1, NA))),
                 "Missing values are not allowed in 'site_cov'.")
    expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                                   repl_cov = list(c = c(1, NA))),
                 "Missing values are not allowed in 'repl_cov'.")
})

