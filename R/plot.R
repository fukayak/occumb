#' @include occumb.R
NULL

setGeneric("plot")

#' Plot method for occumbFit class.
#'
#' @param x An occumbFit object.
#' @param y Not used.
#' @param ... Other arguments passed to the plot method for
#'            \code{jagsUI} class.
#' @export
setMethod("plot", signature(x = "occumbFit"),
    function(x, ...) {
        plot(x@fit, ...)
    }
)

