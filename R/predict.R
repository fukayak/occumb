#' @include occumb.R
NULL

#' @title Predict method for occumbFit class.
#' @description Obtain predictions of parameters related to species occupancy
#'  and detection from an \code{occumbFit} model object.
#' @param object An \code{occumbFit} object.
#' @param newdata An optional \code{occumbData} object with covariates to be
#'  used for prediction. If omitted, the fitted covariates are used.
#' @param parameter The parameter to be predicted.
#' @param scale The scale on which the prediction is made.
#'  \code{type = "response"} returns the prediction on the original scale of the
#'  parameter. \code{type = "link"} returns the prediction on the link scale of
#'  the parameter.
#' @param type The type of prediction.
#'  \code{type = "quantiles"} returns 50% quantile as the posterior median of
#'  the prediction in addition to 2.5 and 97.5% quantiles as the lower and
#'  upper limits of the 95% credible interval of the prediction.
#'  \code{type = "mean"} returns the posterior mean of the prediction.
#'  \code{type = "samples"} returns the posterior samples of the prediction.
#' @param output_dataframe If \code{TRUE}, results are returned in data frame
#'  format.
#' @return
#'  Predictions are obtained as a matrix or array that can have dimensions
#'  corresponding to statistics (or samples), species, sites, and replicates.
#'  The \code{dimension} and \code{label} attributes are added to the output
#'  object to inform these dimensions.
#'  If the sequence read count data \code{y} has species, site, or replicate
#'  names appended as the \code{dimnames} attribute (see Details in
#'  \code{\link{occumbData}()}), they will be copied into the \code{label}
#'  attribute of the returned object.
#'
#'  When \code{output_dataframe = TRUE}, the results are returned in data
#'  frame format where the attributes obtained when
#'  \code{output_dataframe = FALSE} are incorporated into the table.
#' @details
#'  Applying \code{predict()} to an \code{occumbFit} object generates predictions
#'  for the specified parameter (\code{phi}, \code{theta}, or \code{psi}) based
#'  on the estimated effects and the given covariates. It is important to
#'  recognize that the predictions are specific to the individual species being
#'  modeled since they depend on the estimated species-specific effects (i.e.,
#'  \code{alpha}, \code{beta}, and \code{gamma}; see
#'  \href{https://fukayak.github.io/occumb/articles/model_specification.html}{the package vignette} for details).
#'  When providing \code{newdata}, it must thus be assumed that the set of
#'  species contained in \code{newdata} is the same as that of the data being
#'  fitted.
#' @export
setMethod("predict", signature(object = "occumbFit"),
  function(object,
           newdata = NULL,
           parameter = c("phi", "theta", "psi"),
           scale = c("response", "link"),
           type = c("quantiles", "mean", "samples"),
           output_dataframe = FALSE) {

    parameter <- match.arg(parameter)
    scale     <- match.arg(scale)
    type      <- match.arg(type)

    result_array_without_attributes <- get_predict(
      object,
      newdata,
      parameter,
      scale,
      type
    )
    result_array <- add_attributes_predict(
      result_array_without_attributes,
      object,
      newdata,
      parameter,
      scale,
      type
    )

    if (output_dataframe) {
      result_df <- array_to_df_predict(result_array, parameter, scale, type)
      return(result_df)
    } else {
      return(result_array)
    }
  }
)

get_predict <- function(object, newdata, parameter, scale, type) {

  set_inv_link <- function(parameter) {
    switch(parameter,
           phi   = exp,
           theta = stats::plogis,
           psi   = stats::plogis)
  }

  data           <- set_data_predict(newdata, object)
  post_pred_link <- get_post_pred_link(object, data, parameter)
  inv_link       <- set_inv_link(parameter)

  if (type == "quantiles") {
    result <- apply(post_pred_link,
                    2:length(dim(post_pred_link)),
                    stats::quantile, probs = c(0.5, 0.025, 0.975))

    if (scale == "response") {
      result <- inv_link(result)
    }
  }

  if (type == "mean") {
    if (scale == "response") {
      result <- apply(inv_link(post_pred_link),
                      2:length(dim(post_pred_link)),
                      mean)
    }
    if (scale == "link") {
      result <- apply(post_pred_link,
                      2:length(dim(post_pred_link)),
                      mean)
    }

    if (is.null(dim(result))) {
      result <- array(result, c(1, length(result)))
    } else {
      result <- array(result, c(1, dim(result)))
    }
  }

  if (type == "samples") {
    result <- post_pred_link

    if (scale == "response") {
      result <- inv_link(result)
    }
  }

  return(result)
}

set_data_predict <- function(newdata, object) {
  if (missing(newdata) | is.null(newdata)) {
    return(object@data)
  }

  check_newdata(newdata, object)
  data <- newdata
  if (!identical(dimnames(newdata@y)[[1]], dimnames(object@data@y)[[1]])) {
    dimnames(data@y)[[1]] <- dimnames(object@data@y)[[1]]
  }
  return(data)
}

check_newdata <- function(newdata, object) {

  ## Stop if the number of species does not match
  if (!identical(dim(newdata@y)[[1]], dim(object@data@y)[[1]])) {
    stop(sprintf("The number of species in 'newdata' (%s) differs from that in the fitted data (%s).", dim(newdata@y)[[1]], dim(object@data@y)[[1]]))
  }

  ## Stop if covariate names and their order do not match
  name_cov <- name_newcov <- list()
  name_cov$spec_cov <-
    if (is.null(names(object@data@spec_cov))) {
      "(None)"
    } else {
      names(object@data@spec_cov)
    }
  name_newcov$spec_cov <-
    if (is.null(names(newdata@spec_cov))) {
      "(None)"
    } else {
      names(newdata@spec_cov)
    }
  name_cov$site_cov <-
    if (is.null(names(object@data@site_cov))) {
      "(None)"
    } else {
      names(object@data@site_cov)
    }
  name_newcov$site_cov <-
    if (is.null(names(newdata@site_cov))) {
      "(None)"
    } else {
      names(newdata@site_cov)
    }
  name_cov$repl_cov <-
    if (is.null(names(object@data@repl_cov))) {
      "(None)"
    } else {
      names(object@data@repl_cov)
    }
  name_newcov$repl_cov <-
    if (is.null(names(newdata@repl_cov))) {
      "(None)"
    } else {
      names(newdata@repl_cov)
    }
  name_mismatch <- c(!identical(name_cov$spec_cov, name_newcov$spec_cov),
                     !identical(name_cov$site_cov, name_newcov$site_cov),
                     !identical(name_cov$repl_cov, name_newcov$repl_cov))
  if (any(name_mismatch)) {
    stop(paste(c("The names of the covariates in 'newdata' and their order must match those in the fitted data.",
                 sprintf("\n  %s: %s (newdata); %s (fitted data)",
                         c("spec_cov", "site_cov", "repl_cov")[name_mismatch],
                         sapply(name_newcov, paste, collapse = ", ")[name_mismatch],
                         sapply(name_cov, paste, collapse = ", ")[name_mismatch])
               ), collapse = ""))
  }

  ## Stop if covariate classes do not match
  class_cov <- unlist(c(sapply(object@data@spec_cov, class),
                        sapply(object@data@site_cov, class),
                        sapply(object@data@repl_cov, function(x) class(c(x)))))
  class_newcov <- unlist(c(sapply(newdata@spec_cov, class),
                           sapply(newdata@site_cov, class),
                           sapply(newdata@repl_cov, function(x) class(c(x)))))
  class_mismatch <- (class_cov != class_newcov)
  if (any(class_mismatch)) {
    stop(paste(c("The covariate classes in 'newdata' must match those in the fitted data.",
                 sprintf("\n  %s: %s (newdata), %s (fitted data)",
                         names(class_cov[class_mismatch]),
                         class_newcov[class_mismatch],
                         class_cov[class_mismatch])), collapse = ""))
  }

  ## Stop if a discrete covariate in newdata contains a new level
  level_cov <- level_newcov <- list()
  for (i in seq_along(newdata@spec_cov)) {
    if (class(object@data@spec_cov[[i]]) %in% c("character", "factor")) {
      eval(parse(text = sprintf("level_cov$%s <- unique(object@data@spec_cov[[i]])",
                                names(object@data@spec_cov)[i])))
      eval(parse(text = sprintf("level_newcov$%s <- unique(newdata@spec_cov[[i]])",
                                names(newdata@spec_cov)[i])))
    }
  }
  for (i in seq_along(newdata@site_cov)) {
    if (inherits(object@data@site_cov[[i]], c("factor", "character"))) {
      eval(parse(text = sprintf("level_cov$%s <- unique(object@data@site_cov[[i]])",
                                names(object@data@site_cov)[i])))
      eval(parse(text = sprintf("level_newcov$%s <- unique(newdata@site_cov[[i]])",
                                names(newdata@site_cov)[i])))
    }
  }
  for (i in seq_along(newdata@repl_cov)) {
    if (inherits(c(object@data@repl_cov[[i]]), "character")) {
      eval(parse(text = sprintf("level_cov$%s <- unique(object@data@repl_cov[[i]])",
                                names(object@data@repl_cov)[i])))
      eval(parse(text = sprintf("level_newcov$%s <- unique(newdata@repl_cov[[i]])",
                                names(newdata@repl_cov)[i])))
    }
  }
  level_mismatch <- vector(length = length(level_cov))
  for (i in seq_along(level_cov)) {
    level_mismatch[i] <- !identical(level_cov[[i]], level_newcov[[i]])
  }
  if (any(level_mismatch)) {
    stop(paste(c("The levels of discrete covariates in 'newdata' must match those in the fitted data.",
                 sprintf("\n  %s: %s (newdata); %s (fitted data)",
                         names(level_cov)[level_mismatch],
                         sapply(level_newcov, paste, collapse = ", ")[level_mismatch],
                         sapply(level_cov, paste, collapse = ", ")[level_mismatch])
               ), collapse = ""))
  }

  ## Stop if an order of factor levels does not match
  level_cov <- level_newcov <- list()
  for (i in seq_along(newdata@spec_cov)) {
    if (inherits(object@data@spec_cov[[i]], "factor")) {
      eval(parse(text = sprintf("level_cov$%s <- levels(object@data@spec_cov[[i]])",
                                names(object@data@spec_cov)[i])))
      eval(parse(text = sprintf("level_newcov$%s <- levels(newdata@spec_cov[[i]])",
                                names(newdata@spec_cov)[i])))
    }
  }
  for (i in seq_along(newdata@site_cov)) {
    if (inherits(object@data@site_cov[[i]], "factor")) {
      eval(parse(text = sprintf("level_cov$%s <- levels(object@data@site_cov[[i]])",
                                names(object@data@site_cov)[i])))
      eval(parse(text = sprintf("level_newcov$%s <- levels(newdata@site_cov[[i]])",
                                names(newdata@site_cov)[i])))
    }
  }
  level_mismatch <- vector(length = length(level_cov))
  for (i in seq_along(level_cov)) {
    level_mismatch[i] <- !identical(level_cov[[i]], level_newcov[[i]])
  }
  if (any(level_mismatch)) {
    stop(paste(c("The levels of discrete covariates in 'newdata' must match those in the fitted data.",
                 sprintf("\n  %s: %s (newdata); %s (fitted data)",
                         names(level_cov)[level_mismatch],
                         sapply(level_newcov, paste, collapse = ", ")[level_mismatch],
                         sapply(level_cov, paste, collapse = ", ")[level_mismatch])
               ), collapse = ""))
  }

  # Warn if the list of species names does not match
  if (!identical(dimnames(newdata@y)[[1]], dimnames(object@data@y)[[1]])) {
    warning("The list of species names in 'newdata' does not match that in the fitted data; the list of species names in the fitted data will be added to the 'label' attribute of the returned object.")
  }
}


get_post_pred_link <- function(object, data, parameter) {

  # Input:
  #   (a matrix, a vector) ... intercept + covariates, or
  #   (a vector, a numeric) ... intercept only or single covariate (for shared_effect)
  # Output:
  #   a vector
  get_post_linpred <- function(post_effect, cov) {

    if (is.null(dim(post_effect)) && (length(cov) == 1)) {
      result <- post_effect * cov
    } else {
      result <- c(post_effect %*% cov)
    }

    return(result)
  }

  formula <-
    switch(parameter,
           phi = formula(object@occumb_args$formula_phi),
           theta = formula(object@occumb_args$formula_theta),
           psi = formula(object@occumb_args$formula_psi))
  formula_shared <-
    switch(parameter,
           phi = formula(object@occumb_args$formula_phi_shared),
           theta = formula(object@occumb_args$formula_theta_shared),
           psi = formula(object@occumb_args$formula_psi_shared))
  effect <-
    switch(parameter,
           phi = "alpha",
           theta = "beta",
           psi = "gamma")
  effect_shared <-
    switch(parameter,
           phi = "alpha_shared",
           theta = "beta_shared",
           psi = "gamma_shared")

  list_cov <- get_covariates(object, parameter)
  has_shared_effect <- !is.null(list_cov$cov_shared)

  I <- dim(data@y)[1]
  J <- dim(data@y)[2]
  K <- dim(data@y)[3]
  N <- object@fit$mcmc.info$n.samples

  post_effect <- get_post_samples(object, effect)
  if (has_shared_effect) {
    post_effect_shared <- get_post_samples(object, effect_shared)
    if (list_cov$type == "i") {
      pred_link <- matrix(nrow = N, ncol = I)
      for (i in seq_len(I)) {
        pred_link[, i] <-
          get_post_linpred(post_effect[, i, ], list_cov$cov) +
          get_post_linpred(post_effect_shared, list_cov$cov_shared[i, ])
      }
    } else if (list_cov$type == "ij") {
      pred_link <- array(dim = c(N, I, J))
      for (i in seq_len(I)) {
        for (j in seq_len(J)) {
          pred_link[, i, j] <-
            get_post_linpred(post_effect[, i, ], list_cov$cov[j, ]) +
            get_post_linpred(post_effect_shared, list_cov$cov_shared[i, j, ])
        }
      }
    } else if (list_cov$type == "ijk") {
      pred_link <- array(dim = c(N, I, J, K))
      for (i in seq_len(I)) {
        for (j in seq_len(J)) {
          for (k in seq_len(K)) {
            pred_link[, i, j, k] <-
              get_post_linpred(post_effect[, i, ], list_cov$cov[j, k, ]) +
              get_post_linpred(post_effect_shared, list_cov$cov_shared[i, j, k, ])
          }
        }
      }
    }
  } else {
    if (list_cov$type == "i") {
      pred_link <- matrix(nrow = N, ncol = I)
      for (i in seq_len(I)) {
        pred_link[, i] <- get_post_linpred(post_effect[, i, ], list_cov$cov)
      }
    } else if (list_cov$type == "ij") {
      pred_link <- array(dim = c(N, I, J))
      for (i in seq_len(I)) {
        for (j in seq_len(J)) {
          pred_link[, i, j] <-
            get_post_linpred(post_effect[, i, ], list_cov$cov[j, ])
        }
      }
    } else if (list_cov$type == "ijk") {
      pred_link <- array(dim = c(N, I, J, K))
      for (i in seq_len(I)) {
        for (j in seq_len(J)) {
          for (k in seq_len(K)) {
            pred_link[, i, j, k] <-
              get_post_linpred(post_effect[, i, ], list_cov$cov[j, k, ])
          }
        }
      }
    }
  }

  return(pred_link)
}


add_attributes_predict <- function(x, object, newdata, parameter, scale, type) {

  label_pred <- function(x, data) {

    get_dimnames <- function(dn) {
      if (!is.null(dn)) {
        return(dn)
      } else {
        return(NULL)
      }
    }

    if (length(dim(x)) == 2 || is.null(dim(x))) {
      out <- list(Species = get_dimnames(dimnames(data@y)[[1]]))
    } else if (length(dim(x)) == 3) {
      out <- list(Species = get_dimnames(dimnames(data@y)[[1]]),
                  Sites = get_dimnames(dimnames(data@y)[[2]]))
    } else if (length(dim(x)) == 4) {
      out <- list(Species = get_dimnames(dimnames(data@y)[[1]]),
                  Sites = get_dimnames(dimnames(data@y)[[2]]),
                  Replicates = get_dimnames(dimnames(data@y)[[3]]))
    }

    return(out)
  }

  data <- set_data_predict(newdata, object)

  attr(x, "parameter") <- parameter
  attr(x, "scale")     <- scale

  dim_pred <-
    if (length(dim(x)) == 2) {
      c("Species")
    } else if (length(dim(x)) == 3) {
      c("Species", "Sites")
    } else if (length(dim(x)) == 4) {
      c("Species", "Sites", "Replicates")
    }

  if (type == "quantiles") {
    attr(x, "dimension") <- c("Statistics", dim_pred)
    attr(x, "label") <- c(list(Statistics = c("50%", "2.5%", "97.5%")),
                          label_pred(x, data))
  } else if (type == "mean") {
    attr(x, "dimension") <- c("Statistics", dim_pred)
    attr(x, "label") <- c(list(Statistics = c("mean")),
                          label_pred(x, data))
  } else if (type == "samples") {
    attr(x, "dimension") <- c("Samples", dim_pred)
    attr(x, "label") <- c(list(Samples = NULL),
                          label_pred(x, data))
  }

  return(x)
}

array_to_df_predict <- function(x_array, parameter, scale, type) {

  convert_to_df <- function(x_array) {
    x_df <- as.data.frame.table(x_array)
    colnames(x_df) <- c(attributes(x_array)$dimension, "Value")
    return(x_df)
  }

  assign_labels <- function(x_df_without_labels, x_array, start) {
    result <- x_df_without_labels

    for (i in start:(ncol(x_df_without_labels) - 1)) {
      if (is.null(attributes(x_array)$label[[i]])) {
        levels(result[, i]) <- seq_along(levels(x_df_without_labels[, i]))
      } else {
        levels(result[, i]) <- attributes(x_array)$label[[i]]
      }
    }

    return(result)
  }

  if (type == "quantiles") {
    x_df <- assign_labels(convert_to_df(x_array), x_array, 2)
  }

  if (type == "mean") {
    x_df_without_labels <- convert_to_df(x_array)
    x_df_without_labels$Statistics <- factor("mean")
    x_df <- assign_labels(x_df_without_labels, x_array, 2)
  }

  if (type == "samples") {
    x_df <- assign_labels(convert_to_df(x_array), x_array, 1)
  }

  return(data.frame(Parameter = factor(parameter), Scale = factor(scale), x_df))
}
