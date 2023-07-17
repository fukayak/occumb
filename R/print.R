#' @include occumb.R
NULL

setGeneric("print")

setMethod("print", signature(x = "occumbFit"),
    function(x, ...) {
        print(x@fit, ...)
    }
)

