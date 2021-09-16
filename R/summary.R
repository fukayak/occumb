#' @include occumb.R
NULL

setGeneric("summary")

#' @export
setMethod("summary", signature(object = "occumbFit"),
    function(object, ...) {
        summary(object@fit, ...)
    }
)

