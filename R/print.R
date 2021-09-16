#' @include occumb.R
NULL

setGeneric("print")

#' @export
setMethod("print", signature(x = "occumbFit"),
    function(x, ...) {
        print(x@fit, ...)
    }
)

