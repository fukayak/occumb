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

    ## y elements are integers.
    if (sum(object@y %% 1 != 0))
        msg <- c(msg,
                 "'y' contains non-integer value(s).")

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
    if (sum(sapply(object@repl_cov, length) != K))
        msg <- c(msg,
                 sprintf("Length of '%s' should match the number of replicates.",
                         names(object@repl_cov)[sapply(object@repl_cov, length) != K]))

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

# Constructor for occumbData
occumbData <- function(y,
                       spec_cov = NULL,
                       site_cov = NULL,
                       repl_cov = NULL) {

    out <- methods::new("occumbData",
                        y = y,
                        spec_cov = spec_cov,
                        site_cov = site_cov,
                        repl_cov = repl_cov)
    return(out)
}

