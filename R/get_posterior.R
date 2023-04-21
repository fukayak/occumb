#' @title Extract posterior samples of parameters from a model-fit object.
#' @description TODO: to be added.
#' @details TODO: to be added.
#' @param fit An \code{occumbFit} object.
#' @param parameter TODO: to be added.
#' @return TODO: to be added.
#' @examples
#' # to be added
#' @export
get_post_samples <- function(
    fit,
    parameter = c("z", "pi", "phi", "theta", "psi",
                  "alpha", "beta", "gamma",
                  "alpha_shared", "beta_shared", "gamma_shared",
                  "Mu", "sigma", "rho")
    ) {

    assert_occumbFit(fit)
    parameter <- match.arg(parameter)
    if (parameter %!in% names(fit@fit$sims.list))
        stop(paste(parameter, "is not included in the fitted model"))

    out <- .get_post_samples(fit, parameter)
    out
}

#' @title Extract posterior summary of parameters from a model-fit object.
#' @description TODO: to be added.
#' @details TODO: to be added.
#' @param fit An \code{occumbFit} object.
#' @param parameter TODO: to be added.
#' @return TODO: to be added.
#' @examples
#' # to be added
#' @export
get_post_summary <- function(
    fit,
    parameter = c("z", "pi", "phi", "theta", "psi",
                  "alpha", "beta", "gamma",
                  "alpha_shared", "beta_shared", "gamma_shared",
                  "Mu", "sigma", "rho")
    ) {

    assert_occumbFit(fit)
    parameter <- match.arg(parameter)
    if (parameter %!in% names(fit@fit$sims.list))
        stop(paste(parameter, "is not included in the fitted model"))

    out <- .get_post_summary(fit, parameter)
    out
}

.get_post_samples <- function(fit, parameter) {

    add_attributes1 <- function(sims_list, type = c("i", "ij", "ijk")) {
        if (type == "i") {
            attr(sims_list, "dimension") <- c("Sample", "Species")
            if (!is.null(dimnames(fit@data@y)[[1]]))
                attr(sims_list, "label") <- list(
                    Sample  = NULL,
                    Species = dimnames(fit@data@y)[[1]]
                )
        }

        if (type == "ij") {
            attr(sims_list, "dimension") <- c("Sample", "Species", "Site")
            if (!is.null(dimnames(fit@data@y)[[1]]) |
                !is.null(dimnames(fit@data@y)[[2]]))
                attr(sims_list, "label") <- list(
                    Sample  = NULL,
                    Species = dimnames(fit@data@y)[[1]],
                    Site    = dimnames(fit@data@y)[[2]]
                )
        }

        if (type == "ijk") {
            attr(sims_list, "dimension") <- c("Sample", "Species", "Site", "Replicate")
            if (!is.null(dimnames(fit@data@y)[[1]]) |
                !is.null(dimnames(fit@data@y)[[2]]) |
                !is.null(dimnames(fit@data@y)[[3]]))
                attr(sims_list, "label") <- list(
                    Sample    = NULL,
                    Species   = dimnames(fit@data@y)[[1]],
                    Site      = dimnames(fit@data@y)[[2]],
                    Replicate = dimnames(fit@data@y)[[3]]
                )
        }

        return(sims_list)
    }

    add_attributes2 <- function(sims_list, covariate) {
        attr(sims_list, "dimension") <- c("Sample", "Species", "Effects")

        if (identical(covariate, 1)) {
            effect_name <- "(Intercept)"
        } else {
            effect_name <- dimnames(covariate)[[length(dim(covariate))]]
        }

        if (!is.null(dimnames(fit@data@y)[[1]]))
            attr(sims_list, "label") <- list(
                Sample  = NULL,
                Species = dimnames(fit@data@y)[[1]],
                Effects = effect_name
            )

        return(sims_list)
    }

    add_attributes3 <- function(sims_list, covariate) {
        if (is.null(covariate)) {
            invisible()
        } else {
            attr(sims_list, "dimension") <- c("Sample", "Effects")
            attr(sims_list, "label") <- list(
                Sample  = NULL,
                Effects = dimnames(covariate)[[length(dim(covariate))]]
            )
        }

        return(sims_list)
    }

    add_attributes4 <- function(sims_list, margs, is_rho = FALSE) {
        if (identical(margs$cov_phi, 1)) {
            effect_name_phi <- "phi | (Intercept)"
        } else {
            effect_name_phi <- paste(
                "phi |",
                dimnames(margs$cov_phi)[[length(dim(margs$cov_phi))]]
            )
        }

        if (identical(margs$cov_theta, 1)) {
            effect_name_theta <- "theta | (Intercept)"
        } else {
            effect_name_theta <- paste(
                "theta |",
                dimnames(margs$cov_theta)[[length(dim(margs$cov_theta))]]
            )
        }

        if (identical(margs$cov_psi, 1)) {
            effect_name_psi <- "psi | (Intercept)"
        } else {
            effect_name_psi <- paste(
                "psi |",
                dimnames(margs$cov_psi)[[length(dim(margs$cov_psi))]]
            )
        }

        if (is_rho) {
            attr(sims_list, "dimension") <- c("Sample", "Effects 1", "Effects 2")
            attr(sims_list, "label") <- list(
                Sample  = NULL,
                Effects1 = c(effect_name_phi, effect_name_theta, effect_name_psi),
                Effects2 = c(effect_name_phi, effect_name_theta, effect_name_psi)
            )
        } else {
            attr(sims_list, "dimension") <- c("Sample", "Effects")
            attr(sims_list, "label") <- list(
                Sample  = NULL,
                Effects = c(effect_name_phi, effect_name_theta, effect_name_psi)
            )
        }

        return(sims_list)
    }

    out <- eval(parse(text = paste0("fit@fit$sims.list$", parameter)))
    margs <- set_modargs(fit@occumb_args$formula_phi,
                         fit@occumb_args$formula_theta,
                         fit@occumb_args$formula_psi,
                         fit@occumb_args$formula_phi_shared,
                         fit@occumb_args$formula_theta_shared,
                         fit@occumb_args$formula_psi_shared,
                         fit@data)

    if (parameter == "z")
        out <- add_attributes1(out, "ij")

    if (parameter == "pi")
        out <- add_attributes1(out, "ijk")

    if (parameter == "phi")
        out <- add_attributes1(out, margs$phi)

    if (parameter == "theta")
        out <- add_attributes1(out, margs$theta)

    if (parameter == "psi")
        out <- add_attributes1(out, margs$psi)

    if (parameter == "alpha")
        out <- add_attributes2(out, margs$cov_phi)

    if (parameter == "beta")
        out <- add_attributes2(out, margs$cov_theta)

    if (parameter == "gamma")
        out <- add_attributes2(out, margs$cov_psi)

    if (parameter == "alpha_shared")
        out <- add_attributes3(out, margs$cov_phi_shared)

    if (parameter == "beta_shared")
        out <- add_attributes3(out, margs$cov_theta_shared)

    if (parameter == "gamma_shared")
        out <- add_attributes3(out, margs$cov_gamma_shared)

    if (parameter == "Mu")
        out <- add_attributes4(out, margs)

    if (parameter == "sigma")
        out <- add_attributes4(out, margs)

    if (parameter == "rho")
        out <- add_attributes4(out, margs, TRUE)

    out
}

.get_post_summary <- function(fit, parameter) {

    rows_extract <- function(fit, parameter) {
        pattern <- paste0(parameter, "\\[")
        return(grep(pattern, rownames(fit@fit$summary)))
    }

    add_attributes1 <- function(tab_summary, type = c("i", "ij", "ijk")) {
        if (type == "i") {
            attr(tab_summary, "dimension") <- c("Species")
            if (!is.null(dimnames(fit@data@y)[[1]]))
                attr(tab_summary, "label") <- list(
                    Species = dimnames(fit@data@y)[[1]]
                )
        }

        if (type == "ij") {
            attr(tab_summary, "dimension") <- c("Species", "Site")
            if (!is.null(dimnames(fit@data@y)[[1]]) |
                !is.null(dimnames(fit@data@y)[[2]]))
                attr(tab_summary, "label") <- list(
                    Species = dimnames(fit@data@y)[[1]],
                    Site    = dimnames(fit@data@y)[[2]]
                )
        }

        if (type == "ijk") {
            attr(tab_summary, "dimension") <- c("Species", "Site", "Replicate")
            if (!is.null(dimnames(fit@data@y)[[1]]) |
                !is.null(dimnames(fit@data@y)[[2]]) |
                !is.null(dimnames(fit@data@y)[[3]]))
                attr(tab_summary, "label") <- list(
                    Species   = dimnames(fit@data@y)[[1]],
                    Site      = dimnames(fit@data@y)[[2]],
                    Replicate = dimnames(fit@data@y)[[3]]
                )
        }

        return(tab_summary)
    }

    out <- fit@fit$summary[rows_extract(fit, parameter), ]

    if (parameter == "z")
        out <- add_attributes1(out, "ij")

    if (parameter == "pi")
        out <- add_attributes1(out, "ijk")

    out

}

