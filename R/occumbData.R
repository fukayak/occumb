#' @include classes.R
NULL

# Validation for occumbData
validate_occumbData <- function(object) {
    msg <- NULL

    ## y is a 3D-array.
    if (length(dim(object@y)) != 3)
        msg <- c(msg,
                 "'y' should be a 3D-array.")

    I <- dim(object@y)[1] # Number of species
    J <- dim(object@y)[2] # Number of sites
    K <- dim(object@y)[3] # Number of replicates

    if (sum(is.na(object@y)) > 0) {
        ## No missing values in y.
        msg <- c(msg,
                 "Missing values are not allowed in 'y'.")
    } else if (sum(object@y %% 1 != 0)) {
        ## y elements are integers.
        msg <- c(msg,
                 "'y' contains non-integer value(s).")
    }

    ## No overlap in the covariate names.
    cov_names <- c(names(object@spec_cov),
                   names(object@site_cov),
                   names(object@repl_cov))
    if(sum(table(cov_names) > 1))
        msg <- c(msg,
                 sprintf("Duplicated covariate names are not allowed: '%s'",
                         names(table(cov_names))[table(cov_names) > 1]))

    ## Appropriate covariate dimensions.
    if (sum(sapply(object@spec_cov, length) != I))
        msg <- c(msg,
                 sprintf("Length of '%s' should match the number of species.",
                         names(object@spec_cov)[sapply(object@spec_cov, length) != I]))
    if (sum(sapply(object@site_cov, length) != J))
        msg <- c(msg,
                 sprintf("Length of '%s' should match the number of sites.",
                         names(object@site_cov)[sapply(object@site_cov, length) != J]))

    check_repl_cov <- vector(length = length(object@repl_cov))
    for (i in seq_along(object@repl_cov)) {
        if (is.matrix(object@repl_cov[[i]]))
            check_repl_cov[i] <- !identical(dim(object@repl_cov[[i]]), c(J, K))
        else
            check_repl_cov[i] <- TRUE
    }
    if (sum(check_repl_cov))
        msg <- c(msg,
                 sprintf("'%s' should have J rows and K columns.",
                         names(object@repl_cov)[check_repl_cov]))

    ## No missing values in covariates.
    if (sum(is.na(unlist(object@spec_cov))) > 0)
        msg <- c(msg,
                 "Missing values are not allowed in 'spec_cov'.")
    if (sum(is.na(unlist(object@site_cov))) > 0)
        msg <- c(msg,
                 "Missing values are not allowed in 'site_cov'.")
    if (sum(is.na(unlist(object@repl_cov))) > 0)
        msg <- c(msg,
                 "Missing values are not allowed in 'repl_cov'.")

    ifelse(is.null(msg), TRUE, msg)
}

# Data format class for occumb
setClass("occumbData",
         slots = c(y = "array",
                   spec_cov = "optional_list",
                   site_cov = "optional_list",
                   repl_cov = "optional_list"),
         validity = validate_occumbData)

#' Constructor for occumbData data class.
#' 
#' \code{occumbData} creates a data list compatible with the model-fitting
#' function \code{occumb}.
#' 
#' @param y A 3-D array of sequence read counts (integer values).
#'          Dimensions are ordered by species, site, and replicate.
#'          Data for missing replicates must be represented by zero vectors.
#'          NAs are not allowed.
#' @param spec_cov A named list of species covariates.
#'                 Each element must be a vector of numeric, factor, or
#'                 character whose length is equal to dim(y)\[1\] (i.e., the
#'                 number of species). Covariates in character are automatically
#'                 converted to factors. NAs are not allowed.
#' @param site_cov A named list of site covariates.
#'                 Each element must be a vector of numeric, factor, or
#'                 character whose length is equal to dim(y)\[2\] (i.e., the
#'                 number of sites). Covariates in character are automatically
#'                 converted to factors. NAs are not allowed.
#' @param repl_cov A named list of replicate covariates.
#'                 Each element must be a matrix of numeric, factor, or
#'                 character whose dimension is equal to dim(y)\[2:3\] (i.e.,
#'                 the number of sites * number of replicates). Covariates in
#'                 character are automatically converted to factors. NAs are
#'                 not allowed.
#' @section Details:
#'      The element names for spec_cov, site_cov, and repl_cov must all be unique.
#' @return  An S4 object of the `occumbData` class.
#' @examples
#' # Generate a small, random dataset
#' data <- occumbData(
#'     y = array(sample.int(8), dim = rep(2, 3)),
#'     spec_cov = list(cov1 = rnorm(2)),
#'     site_cov = list(cov2 = rnorm(2),
#'                     cov3 = factor(1:2)))
#' @export
occumbData <- function(y,
                       spec_cov = NULL,
                       site_cov = NULL,
                       repl_cov = NULL) {

    # Convert character covariates to factor

    out <- methods::new("occumbData",
                        y = y,
                        spec_cov = spec_cov,
                        site_cov = site_cov,
                        repl_cov = repl_cov)
    return(out)
}

