#' @name get_posterior
#' @title Extract posterior samples or summary of parameters from a model-fit object.
#' @description
#'  \code{get_post_samples()} extracts posterior samples of the specified
#'  parameters from a model-fit object.
#'
#'  \code{get_post_summary()} extracts posterior summary of the specified
#'  parameters from a model-fit object.
#' @details The functions return posterior samples or a summary of one of the
#'  following parameters in the model, stored in the model-fit object
#'  \code{fit}:
#'  \describe{
#'  \item{\code{z}}{Site occupancy status of species.}
#'  \item{\code{pi}}{Multinomial probabilities of species sequence read counts.}
#'  \item{\code{phi}}{Sequence relative dominance of species.}
#'  \item{\code{theta}}{Sequence capture probabilities of species.}
#'  \item{\code{psi}}{Site occupancy probabilities of species.}
#'  \item{\code{alpha}}{Species-specific effects on sequence relative dominance
#'  (\code{phi}).}
#'  \item{\code{beta}}{Species-specific effects on sequence capture
#'  probabilities (\code{theta}).}
#'  \item{\code{gamma}}{Species-specific effects on site occupancy
#'  probabilities (\code{psi}).}
#'  \item{\code{alpha_shared}}{Effects on sequence relative dominance
#'  (\code{phi}) common across species.}
#'  \item{\code{beta_shared}}{Effects on sequence capture probabilities
#'  (\code{theta}) that are common across species.}
#'  \item{\code{gamma_shared}}{Effects on site occupancy probabilities
#'  (\code{psi}) that are common across species.}
#'  \item{\code{Mu}}{Community-level averages of species-specific effects
#'  (\code{alpha}, \code{beta}, \code{gamma}).}
#'  \item{\code{sigma}}{Standard deviations of species-specific effects
#'  (\code{alpha}, \code{beta}, \code{gamma}).}
#'  \item{\code{rho}}{Correlation coefficients of the species-specific effects
#'  (\code{alpha}, \code{beta}, \code{gamma}).}
#'  }
#'  See \href{https://fukayak.github.io/occumb/articles/model_specification.html}{the package vignette}
#'  for details of these parameters.
#'
#'  The parameter may have dimensions corresponding to species, sites,
#'  replicates, and effects (covariates), and when
#'  \code{output_dataframe = FALSE}, the \code{dimension} and \code{label}
#'  attributes are added to the output object to inform these dimensions.
#'  If the sequence read count data \code{y} have species, site, or replicate
#'  names appended as the \code{dimnames} attribute (see Details in
#'  \code{\link{occumbData}()}), they are copied into the \code{label}
#'  attribute of the returned object.
#' @param fit An \code{occumbFit} object.
#' @param parameter A string of parameter name. See Details for possible choices
#'  and corresponding parameters.
#' @param output_dataframe If \code{TRUE}, results are returned in data frame
#'  format.
#' @return
#'  By default, \code{get_post_samples()} returns a vector, matrix, or array of
#'  posterior samples for a selected parameter.
#'
#'  \code{get_post_summary()} returns, by default, a table (matrix) of the
#'  posterior summary of the selected parameters. The elements of the posterior
#'  summary are the same as those obtained with the \code{\link[jagsUI]{jags}()}
#'  function in the \code{jagsUI} package: they include the mean, standard
#'  deviation, percentiles of posterior samples; the \code{Rhat} statistic;
#'  the effective sample size, \code{n.eff}; \code{overlap0}, which checks if 0
#'  falls in the parameter's 95% credible interval; and the proportion of the
#'  posterior with the same sign as the mean, \code{f}.
#'
#'  The \code{dimension} and \code{label} attributes of the output object
#'  provide information regarding the dimensions of the parameter.
#'
#'  When \code{output_dataframe = TRUE}, the results are returned in data
#'  frame format where the attributes obtained when
#'  \code{output_dataframe = FALSE} are incorporated into the table.
#' @examples
#' \donttest{
#' # Generate the smallest random dataset (2 species * 2 sites * 2 reps)
#' I <- 2 # Number of species
#' J <- 2 # Number of sites
#' K <- 2 # Number of replicates
#' y_named <- array(sample.int(I * J * K), dim = c(I, J, K))
#' dimnames(y_named) <- list(c("species 1", "species 2"),
#'                           c("site 1", "site 2"), NULL)
#' data_named <- occumbData(y = y_named)
#'
#' # Fitting a null model
#' fit <- occumb(data = data_named, n.iter = 10100)
#'
#' # Extract posterior samples
#' (post_sample_z <- get_post_samples(fit, "z"))
#' # Look dimensions of the parameter
#' attributes(post_sample_z)
#'
#' # Extract posterior summary
#' (post_summary_z <- get_post_summary(fit, "z"))
#' # Look dimensions of the parameter
#' attributes(post_summary_z)
#' }
NULL

#' @rdname get_posterior
#' @export
get_post_samples <- function(
  fit,
  parameter = c("z", "pi", "phi", "theta", "psi",
                "alpha", "beta", "gamma",
                "alpha_shared", "beta_shared", "gamma_shared",
                "Mu", "sigma", "rho"),
  output_dataframe = FALSE
) {

  # Validate arguments
  parameter <- match.arg(parameter)
  check_args_get_posterior(fit, parameter)

  # Extract
  out <- .get_post_samples(fit, parameter, output_dataframe)
  out
}

#' @rdname get_posterior
#' @export
get_post_summary <- function(
  fit,
  parameter = c("z", "pi", "phi", "theta", "psi",
                "alpha", "beta", "gamma",
                "alpha_shared", "beta_shared", "gamma_shared",
                "Mu", "sigma", "rho"),
  output_dataframe = FALSE
) {

  # Validate arguments
  parameter <- match.arg(parameter)
  check_args_get_posterior(fit, parameter)

  # Extract
  out <- .get_post_summary(fit, parameter, output_dataframe)
  out
}

check_args_get_posterior <- function(fit, parameter) {
  assert_occumbFit(fit)
  if (parameter %!in% names(fit@fit$sims.list))
    stop(paste(parameter, "is not included in the fitted model\n"))
}

.get_post_samples <- function(fit, parameter, output_dataframe) {

  # Function for converting the resulting array to dataframe
  array_to_df <- function(x, parameter) {
    out_df <- as.data.frame.table(x)
    colnames(out_df) <- c(attributes(x)$dimension, "Value")

    levels(out_df[, 1]) <- seq_along(levels(out_df[, 1]))
    for (i in 2:(ncol(out_df) - 1)) {
      if (is.null(attributes(x)$label[[i]])) {
        levels(out_df[, i]) <- seq_along(levels(out_df[, i]))
      } else {
        levels(out_df[, i]) <- attributes(x)$label[[i]]
      }
    }

    data.frame(Parameter = factor(parameter), out_df)
  }

  # Extract samples of the specified parameter
  samples_extracted <- eval(
    parse(text = paste0("fit@fit$sims.list$", parameter))
  )
  if (is.vector(samples_extracted)) {
    samples_extracted <- matrix(samples_extracted, ncol = 1)
  }

  # Add attributes
  out <- add_attributes(samples_extracted, fit, parameter, "samples")

  if (output_dataframe) {
    array_to_df(out, parameter)
  } else {
    out
  }
}


.get_post_summary <- function(fit, parameter, output_dataframe) {

  # Function for converting the resulting array to dataframe
  array_to_df <- function(x, parameter, fit) {

    get_summary_df <- function(x) {
      if (is.null(dim(x))) {
        summary_df <- data.frame(matrix(x, nrow = 1))
        colnames(summary_df) <- names(x)
        rownames(summary_df) <- attributes(x)$parameter
      } else {
        summary_df <- data.frame(x)
        colnames(summary_df) <- colnames(x)
      }
      summary_df
    }

    get_label_df <- function(x, fit, parameter) {
      get_label <- function(i) {
        if (is.null(attributes(x)$label[[i]])) {
          if (attributes(x)$dimension[[i]] == "Species") {
            seq_len(dim(get_data(fit, "y"))[1])
          } else if (attributes(x)$dimension[[i]] == "Site") {
            seq_len(dim(get_data(fit, "y"))[2])
          } else if (attributes(x)$dimension[[i]] == "Replicate") {
            seq_len(dim(get_data(fit, "y"))[3])
          }
        } else {
          attributes(x)$label[[i]]
        }
      }

      if (length(attributes(x)$dimension) == 1) {
        label_df <- expand.grid(get_label(1))
      }
      if (length(attributes(x)$dimension) == 2) {
        if (attributes(x)$parameter == "rho") {
          idx1 <- unlist(
            sapply(seq_along(attributes(x)$label[[1]]), seq)
          )
          idx2 <- unlist(
            sapply(seq_along(attributes(x)$label[[1]]), \(n) rep(n + 1, n))
          )
          label_df <- data.frame(
            attributes(x)$label[[1]][idx1],
            attributes(x)$label[[2]][idx2]
          )
        } else {
          label_df <- expand.grid(get_label(1),
                                  get_label(2))
        }
      }
      if (length(attributes(x)$dimension) == 3) {
        label_df <- expand.grid(get_label(1),
                                get_label(2),
                                get_label(3))
      }
      for (i in seq_len(ncol(label_df))) {
        label_df[[i]] <- factor(label_df[[i]])
      }
      colnames(label_df) <- attributes(x)$dimension

      label_df
    }

    summary_df <- get_summary_df(x)
    label_df <- get_label_df(x, fit, parameter)

    out <- data.frame(Parameter = factor(parameter), label_df, summary_df)
    colnames(out)[colnames(out) == "X2.5."] <- "2.5%"
    colnames(out)[colnames(out) == "X25."] <- "25%"
    colnames(out)[colnames(out) == "X50."] <- "50%"
    colnames(out)[colnames(out) == "X75."] <- "75%"
    colnames(out)[colnames(out) == "X97.5."] <- "97.5%"
    out
  }

  # Identify rows to extract
  rows_extract <- function(fit, parameter) {
    # Get model arguments
    margs <- set_modargs(stats::as.formula(fit@occumb_args$formula_phi),
                         stats::as.formula(fit@occumb_args$formula_theta),
                         stats::as.formula(fit@occumb_args$formula_psi),
                         stats::as.formula(fit@occumb_args$formula_phi_shared),
                         stats::as.formula(fit@occumb_args$formula_theta_shared),
                         stats::as.formula(fit@occumb_args$formula_psi_shared),
                         fit@data)

    pattern <- paste0(parameter, "\\[")
    if (parameter == "alpha_shared") {
      if (margs$M_phi_shared == 1)
        pattern <- paste0(parameter)
    } else if (parameter == "beta_shared") {
      if (margs$M_theta_shared == 1)
        pattern <- paste0(parameter)
    } else if (parameter == "gamma_shared") {
      if (margs$M_psi_shared == 1)
        pattern <- paste0(parameter)
    }

    return(grep(pattern, rownames(fit@fit$summary)))
  }

  # Extract summary of the specified parameter
  summary_extracted <- fit@fit$summary[rows_extract(fit, parameter), ]

  # Add attributes
  out <- add_attributes(summary_extracted, fit, parameter, "summary")

  if (output_dataframe) {
    array_to_df(out, parameter, fit)
  } else {
    out
  }
}

add_attributes <- function(obj, fit, parameter,
                           type = c("samples", "summary")) {

  # Get model arguments
  margs <- set_modargs(stats::as.formula(fit@occumb_args$formula_phi),
                       stats::as.formula(fit@occumb_args$formula_theta),
                       stats::as.formula(fit@occumb_args$formula_psi),
                       stats::as.formula(fit@occumb_args$formula_phi_shared),
                       stats::as.formula(fit@occumb_args$formula_theta_shared),
                       stats::as.formula(fit@occumb_args$formula_psi_shared),
                       fit@data)

  dimnames_y <- dimnames(fit@data@y)

  if (parameter == "z")
    return(add_attributes1(obj, "ij", dimnames_y, type, parameter))

  if (parameter == "pi")
    return(add_attributes1(obj, "ijk", dimnames_y, type, parameter))

  if (parameter == "phi")
    return(add_attributes1(obj, margs$phi, dimnames_y, type, parameter))

  if (parameter == "theta")
    return(add_attributes1(obj, margs$theta, dimnames_y, type, parameter))

  if (parameter == "psi")
    return(add_attributes1(obj, margs$psi, dimnames_y, type, parameter))

  if (parameter == "alpha")
    return(add_attributes2(obj, margs$cov_phi, dimnames_y, type, parameter))

  if (parameter == "beta")
    return(add_attributes2(obj, margs$cov_theta, dimnames_y, type, parameter))

  if (parameter == "gamma")
    return(add_attributes2(obj, margs$cov_psi, dimnames_y, type, parameter))

  if (parameter == "alpha_shared")
    return(add_attributes3(obj, margs$cov_phi_shared, type, parameter))

  if (parameter == "beta_shared")
    return(add_attributes3(obj, margs$cov_theta_shared, type, parameter))

  if (parameter == "gamma_shared")
    return(add_attributes3(obj, margs$cov_psi_shared, type, parameter))

  if (parameter == "Mu")
    return(add_attributes4(obj, fit, FALSE, type, parameter))

  if (parameter == "sigma")
    return(add_attributes4(obj, fit, FALSE, type, parameter))

  if (parameter == "rho")
    return(add_attributes4(obj, fit, TRUE, type, parameter))
}

add_attributes1 <- function(obj, dimension = c("i", "ij", "ijk"),
                            dimnames_y, type, parameter) {

  attr(obj, "parameter") <- parameter

  if (type == "samples") {
    if (dimension == "i") {
      attr(obj, "dimension") <- c("Sample", "Species")
      if (!is.null(dimnames_y[[1]]))
        attr(obj, "label") <- list(
          Sample  = NULL,
          Species = dimnames_y[[1]]
        )
    }

    if (dimension == "ij") {
      attr(obj, "dimension") <- c("Sample", "Species", "Site")
      if (!is.null(dimnames_y[[1]]) || !is.null(dimnames_y[[2]]))
        attr(obj, "label") <- list(
          Sample  = NULL,
          Species = dimnames_y[[1]],
          Site    = dimnames_y[[2]]
        )
    }

    if (dimension == "ijk") {
      attr(obj, "dimension") <- c("Sample", "Species", "Site", "Replicate")
      if (!is.null(dimnames_y[[1]]) ||
            !is.null(dimnames_y[[2]]) ||
            !is.null(dimnames_y[[3]]))
        attr(obj, "label") <- list(
          Sample    = NULL,
          Species   = dimnames_y[[1]],
          Site      = dimnames_y[[2]],
          Replicate = dimnames_y[[3]]
        )
    }
  }

  if (type == "summary") {
    if (dimension == "i") {
      attr(obj, "dimension") <- c("Species")
      if (!is.null(dimnames_y[[1]]))
        attr(obj, "label") <- list(
          Species = dimnames_y[[1]]
        )
    }

    if (dimension == "ij") {
      attr(obj, "dimension") <- c("Species", "Site")
      if (!is.null(dimnames_y[[1]]) || !is.null(dimnames_y[[2]]))
        attr(obj, "label") <- list(
          Species = dimnames_y[[1]],
          Site    = dimnames_y[[2]]
        )
    }

    if (dimension == "ijk") {
      attr(obj, "dimension") <- c("Species", "Site", "Replicate")
      if (!is.null(dimnames_y[[1]]) ||
            !is.null(dimnames_y[[2]]) ||
            !is.null(dimnames_y[[3]]))
        attr(obj, "label") <- list(
          Species   = dimnames_y[[1]],
          Site      = dimnames_y[[2]],
          Replicate = dimnames_y[[3]]
        )
    }
  }

  return(obj)
}

add_attributes2 <- function(obj, covariate, dimnames_y, type, parameter) {

  has_intercept_only <- function(covariate) {
    if (is.array(covariate)) {
      return(dim(covariate)[length(dim(covariate))] == 1)
    } else {
      return(identical(covariate, 1))
    }
  }

  attr(obj, "parameter") <- parameter

  if (type == "samples") {
    attr(obj, "dimension") <- c("Sample", "Species", "Effects")

    if (has_intercept_only(covariate)) {
      effect_name <- "(Intercept)"
    } else {
      effect_name <- dimnames(covariate)[[length(dim(covariate))]]
    }

    if (!is.null(dimnames_y[[1]])) {
      attr(obj, "label") <- list(
        Sample  = NULL,
        Species = dimnames_y[[1]],
        Effects = effect_name
      )
    } else {
      attr(obj, "label") <- list(
        Sample  = NULL,
        Species = NULL,
        Effects = effect_name
      )
    }
  }

  if (type == "summary") {
    attr(obj, "dimension") <- c("Species", "Effects")

    if (has_intercept_only(covariate)) {
      effect_name <- "(Intercept)"
    } else {
      effect_name <- dimnames(covariate)[[length(dim(covariate))]]
    }

    if (!is.null(dimnames_y[[1]])) {
      attr(obj, "label") <- list(
        Species = dimnames_y[[1]],
        Effects = effect_name
      )
    } else {
      attr(obj, "label") <- list(
        Species = NULL,
        Effects = effect_name
      )
    }
  }

  return(obj)
}

add_attributes3 <- function(obj, covariate, type, parameter) {

  attr(obj, "parameter") <- parameter

  if (type == "samples") {
    if (is.null(covariate)) {
      invisible()
    } else {
      attr(obj, "dimension") <- c("Sample", "Effects")
      attr(obj, "label") <- list(
        Sample  = NULL,
        Effects = dimnames(covariate)[[length(dim(covariate))]]
      )
    }
  }

  if (type == "summary") {
    if (is.null(covariate)) {
      invisible()
    } else {
      attr(obj, "dimension") <- c("Effects")
      attr(obj, "label") <- list(
        Effects = dimnames(covariate)[[length(dim(covariate))]]
      )
    }
  }

  return(obj)
}

add_attributes4 <- function(obj, fit, is_rho, type, parameter) {

  # Get model arguments
  margs <- set_modargs(stats::as.formula(fit@occumb_args$formula_phi),
                       stats::as.formula(fit@occumb_args$formula_theta),
                       stats::as.formula(fit@occumb_args$formula_psi),
                       stats::as.formula(fit@occumb_args$formula_phi_shared),
                       stats::as.formula(fit@occumb_args$formula_theta_shared),
                       stats::as.formula(fit@occumb_args$formula_psi_shared),
                       fit@data)

  attr(obj, "parameter") <- parameter

  if (type == "samples") {
    if (identical(fit@occumb_args$formula_phi, "~ 1")) {
      effect_name_phi <- "phi | (Intercept)"
    } else {
      effect_name_phi <- paste(
        "phi |",
        dimnames(margs$cov_phi)[[length(dim(margs$cov_phi))]]
      )
    }

    if (identical(fit@occumb_args$formula_theta, "~ 1")) {
      effect_name_theta <- "theta | (Intercept)"
    } else {
      effect_name_theta <- paste(
        "theta |",
        dimnames(margs$cov_theta)[[length(dim(margs$cov_theta))]]
      )
    }

    if (identical(fit@occumb_args$formula_psi, "~ 1")) {
      effect_name_psi <- "psi | (Intercept)"
    } else {
      effect_name_psi <- paste(
        "psi |",
        dimnames(margs$cov_psi)[[length(dim(margs$cov_psi))]]
      )
    }

    if (is_rho) {
      attr(obj, "dimension") <- c("Sample", "Effects 1", "Effects 2")
      attr(obj, "label") <- list(
        Sample  = NULL,
        Effects1 = c(effect_name_phi, effect_name_theta, effect_name_psi[-length(effect_name_psi)]),
        Effects2 = c(effect_name_phi, effect_name_theta, effect_name_psi)
      )
    } else {
      attr(obj, "dimension") <- c("Sample", "Effects")
      attr(obj, "label") <- list(
        Sample  = NULL,
        Effects = c(effect_name_phi, effect_name_theta, effect_name_psi)
      )
    }
  }

  if (type == "summary") {
    if (identical(fit@occumb_args$formula_phi, "~ 1")) {
      effect_name_phi <- "phi | (Intercept)"
    } else {
      effect_name_phi <- paste(
        "phi |",
        dimnames(margs$cov_phi)[[length(dim(margs$cov_phi))]]
      )
    }

    if (identical(fit@occumb_args$formula_theta, "~ 1")) {
      effect_name_theta <- "theta | (Intercept)"
    } else {
      effect_name_theta <- paste(
        "theta |",
        dimnames(margs$cov_theta)[[length(dim(margs$cov_theta))]]
      )
    }

    if (identical(fit@occumb_args$formula_psi, "~ 1")) {
      effect_name_psi <- "psi | (Intercept)"
    } else {
      effect_name_psi <- paste(
        "psi |",
        dimnames(margs$cov_psi)[[length(dim(margs$cov_psi))]]
      )
    }

    if (is_rho) {
      attr(obj, "dimension") <- c("Effects 1", "Effects 2")
      attr(obj, "label") <- list(
        Effects1 = c(effect_name_phi, effect_name_theta, effect_name_psi[-length(effect_name_psi)]),
        Effects2 = c(effect_name_phi, effect_name_theta, effect_name_psi)
      )
    } else {
      attr(obj, "dimension") <- c("Effects")
      attr(obj, "label") <- list(
        Effects = c(effect_name_phi, effect_name_theta, effect_name_psi)
      )
    }
  }

  return(obj)
}
