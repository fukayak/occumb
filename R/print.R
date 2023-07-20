#' @include occumb.R
NULL

setGeneric("print")

setMethod("print", signature(x = "occumbFit"),
    function(x, ...) {
        print(x@fit, ...)
    }
)

setMethod("print", signature(x = "occumbGof"),
    function(x, ...) {
        cat(crayon::bold('Posterior predictive check for an occumbFit object:\n'))
        cat(' Statistics:',
            ifelse(x@stats == "Freeman_Tukey", "Freeman-Tukey", x@stats),
            '\n')
        cat(' p-value:   ', round(x@p_value, 5), '\n')
        cat(' Discrepancy statistics for observed data:   ',
            round(mean(x@stats_obs), 2), '(mean),',
            round(sd(x@stats_obs), 2), '(sd)',
            '\n')
        cat(' Discrepancy statistics for replicated data: ',
            round(mean(x@stats_rep), 2), '(mean),',
            round(sd(x@stats_rep), 2), '(sd)',
            '\n')
    }
)

