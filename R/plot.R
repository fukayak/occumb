#' @include occumb.R
NULL

setGeneric("plot")

setMethod("plot", signature(x = "occumbFit"),
    function(x, ...) {
        plot(x@fit, ...)
    }
)

