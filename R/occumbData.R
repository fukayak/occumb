#' @include classes.R
NULL

# Validation for occumbData
validate_occumbData <- function(object) {
  msg <- NULL

  ### Tests for sequence read counts
  ## y is not an array of lists.
  if (is.list(object@y)) {
    msg <- c(msg,
             "Elements of 'y' are lists but should be integers")
    return(msg)
  }

  ## y is a 3d-array.
  if (!checkmate::test_array(object@y, d = 3)) {
    msg <- c(msg, "'y' is not a 3D-array")
  }
  if (!length(object@y)) {
    msg <- c(msg, "'y' is an empty array (length(y) = 0)")
  }

  ## No missing values in y.
  if (!checkmate::test_array(object@y, any.missing = FALSE)) {
    msg <- c(msg, "'y' contains missing value(s)")
  }

  ## y elements are integers.
  if (!checkmate::test_array(object@y, mode = "integerish")) {
    if (is.numeric(object@y) && !checkmate::test_numeric(object@y, upper = .Machine$integer.max)) {
      msg <- c(msg, "'y' contains value(s) exceeding maximum integer size")
    } else {
      msg <- c(msg, "'y' contains non-integer value(s)")
    }
  }

  ## y elements are non-negative.
  if (!checkmate::test_integerish(object@y, lower = 0)) {
    msg <- c(msg, "'y' contains negative value(s)")
  }

  ## y elements are all zero.
  if (length(object@y) && checkmate::test_array(object@y, mode = "integerish")) {
    if (!sum(object@y, na.rm = TRUE)) {
      msg <- c(msg, "'y' contains only zero values")
    }
  }

  I <- dim(object@y)[1] # Number of species
  J <- dim(object@y)[2] # Number of sites
  K <- dim(object@y)[3] # Number of replicates

  ### Tests for covariates
  ## Covariates are named list.
  if (!is.null(object@spec_cov) &&
        !checkmate::test_names(names(object@spec_cov))) {
    msg <- c(msg, "'spec_cov' contains unnamed element(s)")
  }
  if (!is.null(object@site_cov) &&
        !checkmate::test_names(names(object@site_cov))) {
    msg <- c(msg, "'site_cov' contains unnamed element(s)")
  }
  if (!is.null(object@repl_cov) &&
        !checkmate::test_names(names(object@repl_cov))) {
    msg <- c(msg, "'repl_cov' contains unnamed element(s)")
  }

  ## No overlap in the covariate names.
  cov_names <- c(names(object@spec_cov),
                 names(object@site_cov),
                 names(object@repl_cov))
  if (sum(table(cov_names) > 1)) {
    msg <- c(msg,
             sprintf("Duplicated covariate names are not allowed: %s",
                     knitr::combine_words(names(table(cov_names))[table(cov_names) > 1],
                                          before = "'", after = "'", and = "")))
  }

  ## Appropriate covariate dimensions.
  if (sum(sapply(object@spec_cov, length) != I)) {
    wrong_spec_cov <- names(object@spec_cov)[sapply(object@spec_cov, length) != I]
    msg <- c(msg,
             sprintf("%s must have a length equal to the number of species",
                     knitr::combine_words(wrong_spec_cov,
                                          before = "'", after = "'")))
  }
  if (sum(sapply(object@site_cov, length) != J)) {
    wrong_site_cov <- names(object@site_cov)[sapply(object@site_cov, length) != J]
    msg <- c(msg,
             sprintf("%s must have a length equal to the number of sites",
                     knitr::combine_words(wrong_site_cov,
                                          before = "'", after = "'")))
  }
  wrong_repl_cov <- vector(length = length(object@repl_cov))
  for (i in seq_along(object@repl_cov)) {
    if (is.matrix(object@repl_cov[[i]])) {
      wrong_repl_cov[i] <- !identical(dim(object@repl_cov[[i]]), c(J, K))
    } else {
      wrong_repl_cov[i] <- TRUE
    }
  }
  if (sum(wrong_repl_cov)) {
    msg <- c(msg,
             sprintf("%s must have a number of rows equal to the number of species and a number of columns equal to the number of sites",
                     knitr::combine_words(names(object@repl_cov)[wrong_repl_cov],
                                          before = "'", after = "'")))
  }

  ## Appropriate covariate classes.
  valid_class <- c("logical", "numeric", "integer", "factor", "character", "NULL")
  cov_class <- c(sapply(object@spec_cov, class),
                 sapply(object@site_cov, class),
                 sapply(object@repl_cov, function(x) class(c(x))))
  if (any(cov_class %!in% valid_class)) {
    msg <- c(
      msg,
      sprintf("%s must be logical, numeric, integer, factor, or character",
              knitr::combine_words(names(cov_class[cov_class %!in% valid_class]),
                                   before = "'", after = "'"))
    )
  }

  ## No missing values in covariates.
  if (!checkmate::test_vector(unlist(object@spec_cov),
                              any.missing = FALSE,
                              null.ok = TRUE)) {
    msg <- c(msg, "'spec_cov' contains missing value(s)")
  }
  if (!checkmate::test_vector(unlist(object@site_cov),
                              any.missing = FALSE,
                              null.ok = TRUE)) {
    msg <- c(msg, "'site_cov' contains missing value(s)")
  }
  if (!checkmate::test_vector(unlist(object@repl_cov),
                              any.missing = FALSE,
                              null.ok = TRUE)) {
    msg <- c(msg, "'repl_cov' contains missing value(s)")
  }

  ## No infinite values in covariates.
  if (checkmate::anyInfinite(unlist(object@spec_cov))) {
    msg <- c(msg, "'spec_cov' contains infinite value(s)")
  }
  if (checkmate::anyInfinite(unlist(object@site_cov))) {
    msg <- c(msg, "'site_cov' contains infinite value(s)")
  }
  if (checkmate::anyInfinite(unlist(object@repl_cov))) {
    msg <- c(msg, "'repl_cov' contains infinite value(s)")
  }

  ifelse(is.null(msg), TRUE, msg)
}

# Data format class for occumb
setClass("occumbData",
         slots = c(y = "array",
                   spec_cov = "list_or_NULL",
                   site_cov = "list_or_NULL",
                   repl_cov = "list_or_NULL"),
         validity = validate_occumbData)

#' @title Constructor for occumbData data class.
#' @description \code{occumbData()} creates a data list compatible with the model fitting
#' function \code{\link{occumb}()}.
#' The element (i.e., covariate) names for \code{spec_cov}, \code{site_cov}, and
#' \code{repl_cov} must all be unique.
#' If \code{y} has a \code{dimnames} attribute, it is retained in the resulting
#' \code{occumbData} object, and can be referenced in subsequent analyses.
#'
#' @param y A 3-D array or a dataframe of sequence read counts
#'          (\code{integer} values). An array's dimensions are ordered by species,
#'          site, and replicate, and may have a \code{dimnames} attribute.
#'          A dataframe's columns are ordered by species, site,
#'          replicate, and sequence read counts.
#'          The data for missing replicates are represented by zero vectors.
#'          \code{NA}s are not allowed.
#' @param spec_cov A named list of species covariates.
#'                 Each covariate can be a vector of continuous (\code{numeric}
#'                 or \code{integer}) or discrete (\code{logical},
#'                 \code{factor}, or \code{character}) variables whose length
#'                 is \code{dim(y)[1]} (i.e., the number of species).
#'                 \code{NA}s are not allowed.
#' @param site_cov A named list of site covariates.
#'                 Each covariate can be a vector of continuous (\code{numeric}
#'                 or \code{integer}) or discrete (\code{logical},
#'                 \code{factor}, or \code{character}) variables whose length
#'                 is \code{dim(y)[1]} (i.e., the number of sites).
#'                 \code{NA}s are not allowed.
#' @param repl_cov A named list of replicate covariates.
#'                 Each covariate can be a matrix of continuous (\code{numeric}
#'                 or \code{integer}) or discrete (\code{logical} or
#'                 \code{character}) variables with dimensions equal to
#'                 \code{dim(y)[2:3]} (i.e., number of sites \eqn{\times}{*}
#'                 number of replicates).
#'                 \code{NA}s are not allowed.
#' @return  An S4 object of the \code{occumbData} class.
#' @examples
#' # Generate the smallest random dataset (2 species * 2 sites * 2 reps)
#' I <- 2 # Number of species
#' J <- 2 # Number of sites
#' K <- 2 # Number of replicates
#' data <- occumbData(
#'     y = array(sample.int(I * J * K), dim = c(I, J, K)),
#'     spec_cov = list(cov1 = rnorm(I)),
#'     site_cov = list(cov2 = rnorm(J), cov3 = factor(1:J)),
#'     repl_cov = list(cov4 = matrix(rnorm(J * K), J, K))
#' )
#'
#' # A case for named y (with species and site names)
#' y_named <- array(sample.int(I * J * K), dim = c(I, J, K))
#' dimnames(y_named) <- list(c("common species", "uncommon species"),
#'                           c("good site", "bad site"), NULL)
#' data_named <- occumbData(
#'     y = y_named,
#'     spec_cov = list(cov1 = rnorm(I)),
#'     site_cov = list(cov2 = rnorm(J), cov3 = factor(1:J)),
#'     repl_cov = list(cov4 = matrix(rnorm(J * K), J, K))
#' )
#' # A real data example
#' data(fish_raw)
#' fish <- occumbData(
#'     y = fish_raw$y,
#'     spec_cov = list(mismatch = fish_raw$mismatch),
#'     site_cov = list(riverbank = fish_raw$riverbank)
#' )
#'
#' # Get an overview of the datasets
#' summary(data)
#' summary(data_named)
#' summary(fish)
#' @export
occumbData <- function(y,
                       spec_cov = NULL,
                       site_cov = NULL,
                       repl_cov = NULL) {

  out <- methods::new("occumbData",
                      y = df_to_array(y),
                      spec_cov = spec_cov,
                      site_cov = site_cov,
                      repl_cov = repl_cov)
  return(out)
}

df_to_array <- function(y) {
  if (!is.data.frame(y)) {
    return(y)
  }

  if (any(duplicated(y[, -4]))) {
    if (any(is.na(y[, -4]))) {
      stop("species/sites/replicates columns contain missing value(s)\n")
    } else {
      stop("duplicate(s) detected in your dataset, only unique values are allowed\n")
    }
  }
  
  species <- unique(y[, 1])
  sites <- unique(y[, 2])
  replicate <- unique(y[, 3])
  I <- length(species)
  J <- length(sites)
  K <- length(replicate)

  if (prod(c(I, J, K)) != nrow(y)) {
    y_expand <- merge(y,
                      expand.grid(species, sites, replicate),
                      all = TRUE)
    y <- replace(y_expand, is.na(y_expand), 0)
  }

  out <- array(NA, dim = c(I, J, K))
  dimnames(out) <- list(species, sites, replicate)

  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      for (k in seq_len(K)) {
        out[i, j, k] <- x[x[, 1] == species[i] && x[, 2] == sites[j] && x[, 3] == replicate[k], 4]
      }
    }
  }

  out
}
