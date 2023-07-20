#' @include occumb.R
NULL

#setGeneric("show")

setMethod("show", "occumbFit",
    function(object) {
        print(object@fit)
    }
)

setMethod("show", "occumbGof",
    function(object) {
        cat(crayon::bold('Posterior predictive check for an occumbFit object:\n'))
        cat(' Statistics:',
            ifelse(object@stats == "Freeman_Tukey", "Freeman-Tukey", object@stats),
            '\n')
        cat(' p-value:   ', round(object@p_value, 5), '\n')
        cat(' Discrepancy statistics for observed data:  ',
            round(mean(object@stats_obs), 2), '(mean),',
            round(stats::sd(object@stats_obs), 2), '(sd)',
            '\n')
        cat(' Discrepancy statistics for replicated data:',
            round(mean(object@stats_rep), 2), '(mean),',
            round(stats::sd(object@stats_rep), 2), '(sd)',
            '\n')
    }
)

