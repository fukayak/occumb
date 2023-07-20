#' @include occumb.R
NULL

setGeneric("plot")

setMethod("plot", signature(x = "occumbFit"),
    function(x, ...) {
        plot(x@fit, ...)
    }
)

setMethod("plot", signature(x = "occumbGof"),
    function(x, ...) {
        plot_gof(x@stats_obs, x@stats_rep, x@stats)
    }
)

