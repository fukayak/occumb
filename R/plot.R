#' @include occumb.R
NULL

#' @title Plot method for occumbFit class.
#' @description Applies \href{https://cran.r-project.org/package=jagsUI}{jagsUI}'s
#'  plot method to an \code{occumbFit} object to draw trace plots
#'  and density plots of MCMC samples of model parameters.
#' @param x An \code{occumbFit} object.
#' @param y \code{NULL}
#' @param ... Additional arguments passed to the plot method for
#'  \href{https://cran.r-project.org/package=jagsUI}{jagsUI} object.
#' @export
setMethod("plot", "occumbFit",
    function(x, y = NULL, ...) {
        plot(x@fit, ...)
    }
)

#' @title Plot method for occumbGof class.
#' @description Draws a scatter plot of fit statistics.
#' @param x An \code{occumbGof} object.
#' @param y \code{NULL}
#' @param ... Additional arguments passed to the default \code{plot} method.
#' @export
setMethod("plot", "occumbGof",
    function(x, y = NULL, ...) {
        plot_gof(x@stats_obs, x@stats_rep, x@stats, ...)
    }
)

