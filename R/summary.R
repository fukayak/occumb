#' @include occumb.R
NULL

#' @title Summary method for occumbFit class.
#' @description Summarizes model fitting result stored in an \code{occumbFit} object.
#' @param object An \code{occumbFit} object.
#' @return  Returns \code{NULL} invisibly.
#' @export
# This code was created on July 17, 2023 by modifying the source code of the
# summary.jagsUI() function in the jagsUI package licensed under the GPLv3.
setMethod("summary", signature(object = "occumbFit"),
    function(object) {
        cat(crayon::bold('Summary for an', class(object), 'object \n\n'))

        N <- apply(object@data@y, c(2, 3), sum)
        reps_per_site <- apply(N, 1, function(x) sum(x > 0))
        cat(crayon::bold('Summary of data:\n'))
        cat(" Number of species, I =", dim(object@data@y)[1], "\n")
        cat(" Number of sites, J =", dim(object@data@y)[2], "\n")
        cat(" Maximum number of replicates per site, K =", dim(object@data@y)[3], "\n")
        cat(" Number of missing observations =", sum(N == 0), "\n")
        cat(" Number of replicates per site:",
            round(mean(reps_per_site), 2), "(average),",
            round(stats::sd(reps_per_site), 2), "(sd)", "\n")
        cat(" Sequencing depth:",
            round(mean(N[N > 0]), 1), "(average),",
            round(stats::sd(N[N > 0]), 1), "(sd)", "\n\n")

        cat(crayon::bold('Model specification:\n'))
        cat(' formula_phi:         ', object@occumb_args$formula_phi, '\n')
        cat(' formula_theta:       ', object@occumb_args$formula_theta, '\n')
        cat(' formula_psi:         ', object@occumb_args$formula_psi, '\n')
        cat(' formula_phi_shared:  ', object@occumb_args$formula_phi_shared, '\n')
        cat(' formula_theta_shared:', object@occumb_args$formula_theta_shared, '\n')
        cat(' formula_psi_shared:  ', object@occumb_args$formula_psi_shared, '\n')
        cat(' prior_prec:          ', object@occumb_args$prior_prec, '\n')
        cat(' prior_ulim:          ', object@occumb_args$prior_ulim, '\n\n')

        cat(crayon::bold('Saved parameters:\n'), object@fit$parameters, '\n\n')

        if (!object@fit$parallel) {
            cat(crayon::bold('MCMC ran for '),
                crayon::bold(object@fit$mcmc.info$elapsed.mins),
                crayon::bold(' minutes at time '),
                crayon::bold(paste(object@fit$run.date)),
                crayon::bold(':\n'), sep = "")
        } else {
            cat(crayon::bold('MCMC ran in parallel for '),
                crayon::bold(object@fit$mcmc.info$elapsed.mins),
                crayon::bold(' minutes at time '),
                crayon::bold(paste(object@fit$run.date)),
                crayon::bold(':\n'), sep = "")
        }
        cat(' For each of', object@fit$mcmc.info$n.chains, 'chains:\n')
        if (all(object@fit$mcmc.info$sufficient.adapt)) {
            cat('  Adaptation:           ',
                mean(object@fit$mcmc.info$n.adapt),
                'iterations (sufficient)\n')
        } else {
            cat('  Adaptation:           ',
                mean(object@fit$mcmc.info$n.adapt),
                'iterations (possibly insufficient)\n')
        }
        cat('  Burn-in:              ',
            object@fit$mcmc.info$n.burnin,
            'iterations\n')
        cat('  Thin rate:            ',
            object@fit$mcmc.info$n.thin,
            'iterations\n')
        cat('  Total chain length:   ',
            object@fit$mcmc.info$n.iter + mean(object@fit$mcmc.info$n.adapt),
            'iterations\n')
        cat('  Posterior sample size:',
            object@fit$mcmc.info$n.samples / object@fit$mcmc.info$n.chains,
            'draws\n\n')

        cat(crayon::bold('Summary of posterior samples: \n'))
        for (i in seq_along(object@fit$parameters)) {
            print_Rhat_neff(object@fit$parameters[i], object)
        }

        if (object@fit$mcmc.info$n.chains == 1) {
            cat('\n')
            cat('**WARNING** Set n.chains > 1 to monitor Rhat and n.eff values.','\n')
        }
    }
)

#' @title Summary method for occumbData class.
#' @description Summarizes dataset stored in an \code{occumbData} object.
#' @param object An \code{occumbData} object.
#' @return  Returns \code{NULL} invisibly.
#' @export
setMethod("summary", signature(object = "occumbData"),
    function(object) {
        N <- apply(object@y, c(2, 3), sum)
        reps_per_site <- apply(N, 1, function(x) sum(x > 0))

        get_list_cov <- function(covariates) {
            if (is.null(names(covariates))) {
                out <- "(None)"
            } else {
                class_cov <- vector(length = length(covariates))
                for (i in seq_along(covariates)) {
                    if (mode(covariates[[i]]) %in% c("logical", "character") |
                        is.factor(covariates[[i]])) {
                        class_cov[i] <- "(categorical)"
                    } else {
                        class_cov[i] <- "(continuous)"
                    }
                }
                out <- paste(names(covariates), class_cov, sep = " ")
            }
            return(out)
        }

        get_list_labels <- function(labels) {
            if (is.null(labels)) {
                return("(None)")
            } else {
                return(labels)
            }
        }

        cat(crayon::bold("Sequence read counts: \n"))
        cat(" Number of species, I =", dim(object@y)[1], "\n")
        cat(" Number of sites, J =", dim(object@y)[2], "\n")
        cat(" Maximum number of replicates per site, K =", dim(object@y)[3], "\n")
        cat(" Number of missing observations =", sum(N == 0), "\n")
        cat(" Number of replicates per site:",
            round(mean(reps_per_site), 2), "(average),",
            round(stats::sd(reps_per_site), 2), "(sd)", "\n")
        cat(" Sequencing depth:",
            round(mean(N[N > 0]), 1), "(average),",
            round(stats::sd(N[N > 0]), 1), "(sd)", "\n\n")

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

print_Rhat_neff <- function(param, object) {

    if (param == "deviance") {
        post_summary <- object@fit$summary[param, ]
    } else {
        post_summary <- get_post_summary(object, param)
    }

    if (is.null(dim(post_summary))) {
        num_per <- 1
        cat(' ', param, ': \n', sep = "")
        if (param != "deviance")
            cat('  Number of parameters:', num_per, '\n')

        if (object@fit$mcmc.info$n.chain > 1) {
            Rhat <- post_summary["Rhat"]
            neff <- post_summary["n.eff"]
            cat('  Rhat: ', round(Rhat, 3), '\n')
            cat('  n.eff:', neff, '\n')
        } else {
            cat('  Rhat:  (not available)', '\n')
            cat('  n.eff: (not available)', '\n')
        }
    } else {
        num_per <- nrow(post_summary)
        cat(' ', param, ': \n', sep = "")
        cat('  Number of parameters:', num_per, '\n')

        if (object@fit$mcmc.info$n.chain > 1) {
            Rhat <- post_summary[, "Rhat"]
            neff <- post_summary[, "n.eff"]
            if (any(is.na(Rhat))) {
                if (sum(is.na(Rhat)) == length(Rhat)) {
                    cat('  Rhat: ', NA, '(min),', NA, '(median),',
                        NA, '(mean),', NA, '(max),',
                        sum(is.na(Rhat)), '(Number of NAs)', '\n')
                } else {
                    cat('  Rhat: ',
                        round(min(Rhat, na.rm = TRUE), 3), '(min),',
                        round(stats::median(Rhat, na.rm = TRUE), 3), '(median),',
                        round(mean(Rhat, na.rm = TRUE), 3), '(mean),',
                        round(max(Rhat, na.rm = TRUE), 3), '(max),',
                        sum(is.na(Rhat)), '(Number of NAs)', '\n')
                }
            } else {
                cat('  Rhat: ',
                    round(min(Rhat), 3), '(min),',
                    round(stats::median(Rhat), 3), '(median),',
                    round(mean(Rhat), 3), '(mean),',
                    round(max(Rhat), 3), '(max)', '\n')
            }
            cat('  n.eff:',
                round(min(neff), 1), '(min),',
                round(stats::median(neff), 1), '(median),',
                round(mean(neff), 1), '(mean),',
                round(max(neff), 1), '(max)', '\n')
        } else {
            cat('  Rhat:  (not available)', '\n')
            cat('  n.eff: (not available)', '\n')
        }
    }
}

