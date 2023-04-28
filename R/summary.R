#' @include occumb.R
NULL

setGeneric("summary")

setMethod("summary", signature(object = "occumbFit"),
    function(object, ...) {
        summary(object@fit, ...)
    }
)

setMethod("summary", signature(object = "occumbData"),
    function(object, ...) {
        N <- apply(object@y, c(2, 3), sum)
        reps_per_site <- apply(N, 1, function(x) sum(x > 0))

        get_list_cov <- function(covariates) {
            if (is.null(names(covariates))) {
                out <- "(None)"
            } else {
                class_cov <- vector(length = length(covariates))
                for (i in seq_along(covariates)) {
                    if (mode(covariates[[i]]) %in% c("logical", "character") |
                        is.factor(covariates[[i]]))
                        class_cov[i] <- "(categorical)"
                    else
                        class_cov[i] <- "(continuous)"
                }
                out <- paste(names(covariates), class_cov, sep = " ")
            }
            return(out)
        }

        get_list_labels <- function(labels) {
            if (is.null(labels))
                return("(None)")
            else
                return(labels)
        }

        cat(crayon::bold("Sequence read conts: \n"))
        cat("Number of species, I =", dim(object@y)[1], "\n")
        cat("Number of sites, J =", dim(object@y)[2], "\n")
        cat("Maximum number of replicates per site, K =", dim(object@y)[3], "\n")
        cat("Number of missing observations =", sum(N == 0), "\n")
        cat("Number of replicates per site:", mean(reps_per_site), "(average),",
            stats::sd(reps_per_site), "(sd)", "\n")
        cat("Sequencing depth:", mean(N[N > 0]), "(average),",
            stats::sd(N[N > 0]), "(sd)", "\n\n")

        cat(crayon::bold("Species covariates: \n"),
            paste(get_list_cov(object@spec_cov), collapse = ", "),
            "\n")
        cat(crayon::bold("Site covariates: \n"),
            paste(get_list_cov(object@site_cov), collapse = ", "),
            "\n")
        cat(crayon::bold("Replicate covariates: \n"),
            paste(get_list_cov(object@repl_cov), collapse = ", "),
            "\n\n")

        cat(crayon::bold("Labels for species: \n"),
            paste(get_list_labels(dimnames(object@y)[[1]]), collapse = ", "),
            "\n")
        cat(crayon::bold("Labels for sites: \n"),
            paste(get_list_labels(dimnames(object@y)[[2]]), collapse = ", "),
            "\n")
        cat(crayon::bold("Labels for replicates: \n"),
            paste(get_list_labels(dimnames(object@y)[[3]]), collapse = ", "),
            "\n")
    }
)

