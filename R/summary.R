#' @include occumb.R
NULL

setGeneric("summary")

setMethod("summary", signature(object = "occumbFit"),
    function(object, ...) {
        summary(object@fit, ...)
    }
)

