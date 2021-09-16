#' @include occumb.R
NULL

setGeneric("summary")

#' Summary method for occumbFit class.
#'
#' @param object An occumbFit object.
#' @param ... Other arguments passed to the summary method for
#'            \code{jagsUI} class.
#' @export
setMethod("summary", signature(object = "occumbFit"),
    function(object, ...) {
        summary(object@fit, ...)
    }
)

