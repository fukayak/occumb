#' @include occumb.R
NULL

setGeneric("print")

#' Print method for occumbFit class.
#'
#' @param x An occumbFit object.
#' @param ... Other arguments passed to the print method for
#'            \code{jagsUI} class.
#' @export
setMethod("print", signature(x = "occumbFit"),
    function(x, ...) {
        print(x@fit, ...)
    }
)

