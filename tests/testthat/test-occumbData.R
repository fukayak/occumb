### Test data ------------------------------------------------------------------
I <- 2
J <- 2
K <- 2
y <- array(sample.int(I * J * K), dim = c(I, J, K))
dimnames(y) <- list(c("A", "B"),
                    c("a", "b"),
                    c("1", "2"))
df <- as.data.frame.table(y)

### Test for validate_occumbData() ---------------------------------------------
## Tests for sequence read counts

test_that("Mode check for y works", {
  expect_error(new("occumbData", y = array(list(), dim = rep(2, 3))),
               "Elements of 'y' are lists but should be integers")
})

test_that("Check for a dataframe for y works", {
  expect_identical(occumbData(y = y),
                   occumbData(y = df))
  expect_identical(occumbData(y = y),
                   occumbData(y = tibble::tibble(df)))
})

test_that("Dimension check for y works", {
  expect_error(new("occumbData", y = array(1:4, dim = rep(2, 1))),
               "'y' is not a 3D-array")
  expect_error(new("occumbData", y = array(1:4, dim = rep(2, 2))),
               "'y' is not a 3D-array")
  expect_error(new("occumbData", y = array(1:16, dim = rep(2, 4))),
               "'y' is not a 3D-array")
  expect_error(new("occumbData", y = array(numeric(), dim = c(2, 2, 0))),
               "'y' is an empty array")
})

test_that("Check for missing values for y works", {
  expect_error(new("occumbData", y = array(c(NA, 1:7), dim = rep(2, 3))),
               "'y' contains missing value.")
})

test_that("Integer check for y works", {
  expect_error(new("occumbData", y = array(1:8 + 0.1, dim = rep(2, 3))),
               "'y' contains non-integer value")
  expect_error(new("occumbData", y = array(c("1", 1:7), dim = rep(2, 3))),
               "'y' contains non-integer value")
  expect_error(new("occumbData", y = array(c(Inf, 1:7), dim = rep(2, 3))),
               "'y' contains value(s) exceeding maximum integer size",
               fixed = TRUE)
})

test_that("Check for non-negative values for y works", {
  expect_error(new("occumbData", y = array(c(-1, 1:7), dim = rep(2, 3))),
               "'y' contains negative value")
})

test_that("Check for non-zero values for y works", {
  expect_error(new("occumbData", y = array(0, dim = rep(2, 3))),
               "'y' contains only zero values")
})

## Tests for covariates
test_that("Check for covariate names works", {
  unnamed_list1 <- unnamed_list2 <- list(rep(1, 2), rep(1, 2))
  names(unnamed_list2) <- c("a", "")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   spec_cov = unnamed_list1),
               "'spec_cov' contains unnamed element")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   spec_cov = unnamed_list2),
               "'spec_cov' contains unnamed element")
  unnamed_list1 <- unnamed_list2 <- list(rep(1, 3), rep(1, 3))
  names(unnamed_list2) <- c("a", "")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   site_cov = unnamed_list1),
               "'site_cov' contains unnamed element")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   site_cov = unnamed_list2),
               "'site_cov' contains unnamed element")
  unnamed_list1 <- unnamed_list2 <-
    list(matrix(1:2, 1, 2), matrix(1:2, 1, 2))
  names(unnamed_list2) <- c("a", "")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   repl_cov = unnamed_list1),
               "'repl_cov' contains unnamed element")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   repl_cov = unnamed_list2),
               "'repl_cov' contains unnamed element")
})

test_that("Check for covariate name overlap works", {
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   spec_cov = list(a = NULL, b = NULL),
                   site_cov = list(a = NULL)),
               "Duplicated covariate names are not allowed: 'a'")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   spec_cov = list(a = NULL, b = NULL),
                   site_cov = list(a = NULL, b = NULL)),
               "Duplicated covariate names are not allowed: 'a', 'b'")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   spec_cov = list(a = NULL, b = NULL),
                   repl_cov = list(a = NULL)),
               "Duplicated covariate names are not allowed: 'a'")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   spec_cov = list(a = NULL, b = NULL),
                   repl_cov = list(a = NULL, b = NULL)),
               "Duplicated covariate names are not allowed: 'a', 'b'")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   site_cov = list(a = NULL, b = NULL),
                   repl_cov = list(a = NULL)),
               "Duplicated covariate names are not allowed: 'a'")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   site_cov = list(a = NULL, b = NULL),
                   repl_cov = list(a = NULL, b = NULL)),
               "Duplicated covariate names are not allowed: 'a', 'b'")
})

test_that("Dimension check for covariates works", {
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   spec_cov = list(a = rep(1, 1))),
               "'a' must have a length equal to the number of species")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   spec_cov = list(a = rep(1, 1),
                                   b = rep(1, 3))),
               "'a' and 'b' must have a length equal to the number of species")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   site_cov = list(a = rep(1, 1))),
               "'a' must have a length equal to the number of sites")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   site_cov = list(a = rep(1, 1),
                                   b = rep(1, 3))),
               "'a' and 'b' must have a length equal to the number of sites")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   repl_cov = list(a = rep(1, 1))),
               "'a' must have a number of rows equal to the number of species and a number of columns equal to the number of sites")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   repl_cov = list(a = rep(1, 1),
                                   b = rep(1, 3))),
               "'a' and 'b' must have a number of rows equal to the number of species and a number of columns equal to the number of sites")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   repl_cov = list(a = matrix(1:2, 1, 2))),
               "'a' must have a number of rows equal to the number of species and a number of columns equal to the number of sites")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   repl_cov = list(a = matrix(1:2, 1, 2),
                                   b = matrix(1:2, 2, 1))),
               "'a' and 'b' must have a number of rows equal to the number of species and a number of columns equal to the number of sites")
})

test_that("Check for covariate classes works", {
  # complex
  spec_cov <- list(a = rnorm(2), b = rep(0 + 1i, 2))
  site_cov <- list(c = rnorm(2), d = rep(0 + 1i, 2))
  repl_cov <- list(e = matrix(1:4, 2, 2), f = matrix(rep(0 + 1i, 4), 2, 2))
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   spec_cov = spec_cov,
                   site_cov = site_cov,
                   repl_cov = repl_cov),
               "'b', 'd', and 'f' must be logical, numeric, integer, factor, or character")
})

test_that("Check for missing covariate values works", {
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   spec_cov = list(a = c(1, NA))),
               "'spec_cov' contains missing value")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   site_cov = list(b = c(1, NA))),
               "'site_cov' contains missing value")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   repl_cov = list(c = matrix(c(1:3, NA), 2, 2))),
               "'repl_cov' contains missing value")
})

test_that("Check for infinite covariate values works", {
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   spec_cov = list(a = c(1, Inf))),
               "'spec_cov' contains infinite value")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   site_cov = list(b = c(1, Inf))),
               "'site_cov' contains infinite value")
  expect_error(new("occumbData", y = array(1:8, dim = rep(2, 3)),
                   repl_cov = list(c = matrix(c(1:3, Inf), 2, 2))),
               "'repl_cov' contains infinite value")
})

### Test for df_to_array() ---------------------------------------------
test_that("Non-dataframe input is returned itself", {
  expect_identical(y, df_to_array(y))
  expect_identical(as.list(y),
                   df_to_array(as.list(y)))
})

test_that("Check for missing combination works", {
  data_missing <- subset(df, !(Var1 == "A" & Var2 == "a" & Var3 == "1"))
  expect_equal(suppressMessages(df_to_array(data_missing))["A", "a", "1"], 0)
  expect_message(df_to_array(data_missing),
                 "The dataset contained missing obervation(s). Read counts of 0 were assigned to them.",
                 fixed = TRUE)
})

test_that("Check for duplicates works", {
  data_duplicated <- df
  data_duplicated[2, ] <- data_duplicated[1, ]
  expect_error(
    expect_output(df_to_array(data_duplicated),
                  data_duplicated[1:2, ]),
    "The dataset contains duplicate observation(s) listed above. Ensure that the dataset has only unique observations.\n",
    fixed = TRUE
  )
  data_duplicated[2, 4] <- 9999
  expect_error(
    expect_output(df_to_array(data_duplicated),
                  data_duplicated[1:2, ]),
    "The dataset contains duplicate observation(s) listed above. Ensure that the dataset has only unique observations.\n",
    fixed = TRUE
  )
})

test_that("Check for missing values in species/sites/replicates column works", {
  data_NA <- df
  data_NA[1, 3] <- NA
  expect_error(df_to_array(data_NA),
               "NAs are not allowed in the dataset.",
               fixed = TRUE)
})

test_that("df_to_array() processes numeric columns adequately", {
  data_numeric <- df
  data_numeric[, 3] <- as.numeric(data_numeric[, 3])
  expect_identical(y, df_to_array(data_numeric))
})

test_that("df_to_array() processes character columns adeuately", {
  data_character <- df
  data_character[, 1] <- as.character(data_character[, 1])
  expect_identical(y, df_to_array(data_character))
})
