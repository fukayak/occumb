#' @include occumb.R
NULL

setGeneric("plot")

#' @export
setMethod("plot", signature(x = "occumbFit"),
    function(x, ...) {
        plot(x@fit, ...)
    }
)

