#' @include classes.R
NULL

# Validation for occumbData
validate_occumbData <- function(object) {
    msg <- NULL

    if (is.list(object@y)) {
        ## y is not an array of lists.
        msg <- c(msg,
                 "'y' should be a 3D-array of integers, not lists.")
        return(msg)
    } else if (length(dim(object@y)) != 3) {
        ## y is a 3D-array.
        msg <- c(msg,
                 "'y' should be a 3D-array.")
    }

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
    wrong_repl_cov <- vector(length = length(object@repl_cov))
    for (i in seq_along(object@repl_cov)) {
        if (is.matrix(object@repl_cov[[i]])) {
            wrong_repl_cov[i] <- !identical(dim(object@repl_cov[[i]]), c(J, K))
        } else {
            wrong_repl_cov[i] <- TRUE
        }
    }
    if (sum(wrong_repl_cov))
        msg <- c(msg,
                 sprintf("'%s' should be a matrix with J rows and K columns.",
                         names(object@repl_cov)[wrong_repl_cov]))

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
                   spec_cov = "list_or_NULL",
                   site_cov = "list_or_NULL",
                   repl_cov = "list_or_NULL"),
         validity = validate_occumbData)

#' Constructor for occumbData data class.
#' 
#' \code{occumbData()} creates a data list compatible with the model-fitting
#' function \code{\link{occumb}()}.
#' 
#' The element names for \code{spec_cov}, \code{site_cov}, and \code{repl_cov}
#' must all be unique.
#' If \code{y} has a \code{dimnames} attribute, it is retained in the resulting
#' \code{occumbData} object and can be referenced in subsequent analyses.
#'
#' @param y A 3-D array of sequence read counts (integer values) that may have
#'          a \code{dimnames} attribute.
#'          Dimensions are ordered by species, site, and replicate.
#'          Data for missing replicates must be represented by zero vectors.
#'          NAs are not allowed.
#' @param spec_cov A named list of species covariates.
#'                 Each element must be a vector of numeric, factor, or
#'                 character whose length is equal to \code{dim(y)[1]} (i.e.,
#'                 the number of species). Characters are automatically 
#'                 converted to factors. NAs are not allowed.
#' @param site_cov A named list of site covariates.
#'                 Each element must be a vector of numeric, factor, or
#'                 character whose length is equal to \code{dim(y)[2]} (i.e.,
#'                 the number of sites). Characters are automatically converted
#'                 to factors. NAs are not allowed.
#' @param repl_cov A named list of replicate covariates.
#'                 Each element must be a matrix of numeric, factor, or
#'                 character whose dimension is equal to \code{dim(y)[2:3]}
#'                 (i.e., the number of sites \eqn{\times}{*} number of
#'                 replicates). Characters are automatically converted to
#'                 factors. NAs are not allowed.
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
#' 
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

    # The mode of covariates is numeric, factor, or character.
    check_covariate_mode(spec_cov, site_cov, repl_cov)

    out <- methods::new("occumbData",
                        y = y,
                        spec_cov = spec_cov,
                        site_cov = site_cov,
                        repl_cov = repl_cov)
    return(out)
}

# Check for the mode of each covariate
check_covariate_mode <- function(spec_cov, site_cov, repl_cov) {
    wrong_spec_cov_mode <- vector(length = length(spec_cov))
    for (i in seq_along(spec_cov)) {
        if (!(mode(spec_cov[[i]]) == "numeric" |
              mode(spec_cov[[i]]) == "factor" | 
              mode(spec_cov[[i]]) == "character"))
            wrong_spec_cov_mode[i] <- TRUE
    }
    wrong_site_cov_mode <- vector(length = length(site_cov))
    for (i in seq_along(site_cov)) {
        if (!(mode(site_cov[[i]]) == "numeric" |
              mode(site_cov[[i]]) == "factor" | 
              mode(site_cov[[i]]) == "character"))
            wrong_site_cov_mode[i] <- TRUE
    }
    wrong_repl_cov_mode <- vector(length = length(repl_cov))
    for (i in seq_along(repl_cov)) {
        if (!(mode(repl_cov[[i]]) == "numeric" |
              mode(repl_cov[[i]]) == "factor" | 
              mode(repl_cov[[i]]) == "character"))
            wrong_repl_cov_mode[i] <- TRUE
    }

    if (sum(c(wrong_spec_cov_mode, wrong_site_cov_mode, wrong_repl_cov_mode))) {
        var <- c(names(spec_cov)[wrong_spec_cov_mode],
                 names(site_cov)[wrong_site_cov_mode],
                 names(repl_cov)[wrong_repl_cov_mode])
        mod <- NULL
        for (i in seq_len(sum(wrong_spec_cov_mode)))
            mod <- c(mod, mode(spec_cov[[which(wrong_spec_cov_mode)[i]]]))
        for (i in seq_len(sum(wrong_site_cov_mode)))
            mod <- c(mod, mode(site_cov[[which(wrong_site_cov_mode)[i]]]))
        for (i in seq_len(sum(wrong_repl_cov_mode)))
            mod <- c(mod, mode(repl_cov[[which(wrong_repl_cov_mode)[i]]]))

        stop(message = sprintf("Unacceptable mode: the following covariates must be numeric, factor, or character. \n %s",
                               paste(sprintf("%s: %s", var, mod), collapse = "; ")))
    }
}

