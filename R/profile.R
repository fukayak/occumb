#' @title Expected utility for local species diversity assessments.
#'
#' @description `eval_util_L()` evaluates the expected utility for the local
#'  species diversity assessment using Monte Carlo integration.
#' @details
#' The utility of the local species diversity assessment can be defined, for a
#'  given set of sites, as the expected number of detected species per site
#'  (Fukaya et al. 2022). `eval_util_L()` evaluates this utility for the set of
#'  sites included in the `occumbFit` object (specifically, the sites for which
#'  the latent occupancy state \eqn{z}{z} is estimated) for the combination of
#'  `K` and `N` values specified in the `conditions` argument.
#'  Such evaluations can be used to balance `K` and `N` to maximize the utility
#'  under a constant budget (possible combinations of `K` and `N` under a
#'  budget level are easily obtained using `x`; see the example below).
#'  It is also possible to examine how the utility varies with different `K`
#'  and `N` values without setting a budget level, which may be useful to know
#'  the satisfactory level of `K` and `N` from a purely technical point of view.
#' 
#' The expected utility is defined as the expected value of the conditional
#'  utility of the form:
#'  \deqn{U(K, N \mid \boldsymbol{r}, \boldsymbol{u}) = \frac{1}{J}\sum_{j = 1}^{J}\sum_{i = 1}^{I}\left\{1 - \prod_{k = 1}^{K}\left(1 - \frac{u_{ijk}r_{ijk}}{\sum_{m = 1}^{I}u_{mjk}r_{mjk}} \right)^N \right\}}{U(K, N | r, u) = (1 / J) * sum_{j, i}((1 - \prod_{k}(1 - (u[i, j, k] * r[i, j, k])/sum(u[, j, k] * r[, j, k])))^N)}
#'  where \eqn{u_{ijk}}{u[i, j, k]} is a latent indicator variable representing
#'  the inclusion of the sequence of species \eqn{i}{i} in replicate \eqn{k}{k}
#'  at site \eqn{j}{j} and \eqn{r_{ijk}}{r[i, j, k]} is a latent variable that
#'  is proportional to the relative frequency of the sequence of species
#'  \eqn{i}{i}, conditional on its presence in replicate \eqn{k}{k} at site
#'  \eqn{j}{j} (Fukaya et al. 2022).
#' Expectations are taken with respect to the posterior distributions of
#'  \eqn{r}{r} and \eqn{u}{u}, which are evaluated numerically by Monte Carlo
#'  integration using MCMC samples in the `occumbFit` object. Higher
#'  approximation accuracy can be obtained by increasing the value of `N_rep`.
#'
#' The expected utility is evaluated assuming that all replicates are
#'  homogeneous in the sense that the model parameters associated with the
#'  detection process are constant across replicates. For this reason, in the
#'  current version, `eval_util_L()` cannot be applied if the supplied
#'  `occumbFit` object contains a model with replicate covariates.
#'
#' Monte Carlo integration is executed in parallel on multiple CPU cores where
#'  the `cores` argument controls the degree of parallelization. By default, all
#'  cores available in the user's environment are used.
#' @param settings A data frame that specifies a set of conditions under which
#'  the utility is evaluated. It must include a column named `K` and `N`, which
#'  specifies the number of replicates per site and the sequencing depth per
#'  replicate, respectively.
#'  `K` and `N` must be numeric vectors greater than 0. When `K` contains a
#'  decimal, the decimal part is discarded and treated as an integer.
#'  Additional columns are ignored but may be included.
#' @param fit An `occumbFit` object containing a posterior sample of the
#'  relevant parameters.
#' @param N_rep Controls the sample size for Monte Carlo integration.
#'   The integral is evaluated using a total of `N_sample * N_rep` random samples,
#'   where `N_sample` is the size of the MCMC sample provided as `fit`.
#' @param cores The number of cores to use for parallelization.
#' @return A data frame with a column named `Utility` in which the estimates of
#'  expected utility are stored. This is obtained by adding the `Utility` column
#'  to the data frame provided in the `settings` argument.
#' @section References:
#'      K. Fukaya, N. I. Kondo, S. S. Matsuzaki and T. Kadoya (2022)
#'      Multispecies site occupancy modelling and study design for spatially
#'      replicated environmental DNA metabarcoding. \emph{Methods in Ecology
#'      and Evolution} \strong{13}:183--193.
#'      \url{https://doi.org/10.1111/2041-210X.13732}
#' @export
eval_util_L <- function(settings,
                        fit,
                        N_rep = 1,
                        cores = parallel::detectCores()) {

    # Validate arguments
    qc_eval_util_L(settings, fit)

    # Extract posterior samples
    z     <- get_post_samples(fit, "z")
    theta <- get_post_samples(fit, "theta")
    phi   <- get_post_samples(fit, "phi")

    # Adapt theta/phi when they are species-specific
    if (length(dim(theta)) == 2) {
        theta <- array(dim = dim(z))
        for (j in seq_len(dim(theta)[3]))
            theta[, , j] <- get_post_samples(fit, "theta")
    }
    if (length(dim(phi)) == 2) {
        phi <- array(dim = dim(z))
        for (j in seq_len(dim(phi)[3]))
            phi[, , j] <- get_post_samples(fit, "phi")
    }

    # Adapt theta/phi when they are replicate-specific
    # ... to be added

    # Calculate expected utility
    result <- rep(NA, nrow(settings))
    for (i in seq_len(nrow(settings))) {
        result[i] <- eutil(z = z, theta = theta, phi = phi,
                           K = settings[i, "K"], N = settings[i, "N"],
                           scale = "local",
                           N_rep = N_rep, cores = cores)
    }

    # Output
    out <- cbind(settings, result)
    colnames(out)[ncol(out)] <- "Utility"
    out
}

# @title Tentative: the main function (regional). 
# @param settings A dataframe of J, K, and N.
# @param fit An occumbFit object.
# @param N_rep Controls the sample size for Monte Carlo simulation.
#   The integral is evaluated using a total of N_sample * N_rep random samples,
#   where N_sample is the size of the MCMC sample given as `fit.`
# @param core The number of cores to use for parallelization.
# @return The expected utility.
eval_util_R <- function(settings,
                        fit,
                        N_rep = 1,
                        cores = parallel::detectCores()) {

    # Validate arguments ... to be added

    result <- rep(NA, nrow(settings))

    # Extract posterior samples
    psi   <- get_post_samples(fit, "psi")
#    theta <- get_post_samples(fit, "theta")
#    phi   <- get_post_samples(fit, "phi")

    # Adapt psi/theta/phi when they are site- or replicate-specific
    # ... to be added

    # Calculate expected utility
    for (i in seq_len(nrow(settings))) {
        # Generate z (dim = N_sample * N_species * N_site)
        z <- array(NA, dim = c(dim(psi)[1], dim(psi)[2], settings$J[i]))
        for (j in seq_len(dim(z)[3]))
            z[, , j] <- rbinom(dim(psi)[1] * dim(psi)[2], 1, psi)

        # Adapt dimension of theta/phi
        theta <- array(dim = dim(z))
        for (j in seq_len(dim(theta)[3]))
            theta[, , j] <- get_post_samples(fit, "theta")
        phi <- array(dim = dim(z))
        for (j in seq_len(dim(phi)[3]))
            phi[, , j] <- get_post_samples(fit, "phi")

        result[i] <- eutil(z = z, theta = theta, phi = phi,
                           K = settings$K[i], N = settings$N[i],
                           scale = "regional",
                           N_rep = N_rep, cores = cores)
    }

    # Output
    out <- cbind(settings, result)
    colnames(out)[ncol(out)] <- "Utility"
    out
}

qc_eval_util_L <- function(settings, fit) {
    # Assert that settings is a data frame and contains the required columns
    checkmate::assert_data_frame(settings)
    if (!checkmate::testSubset("K", names(settings)))
        stop("The 'settings' argument does not contain column 'K'.")
    if (!checkmate::testSubset("N", names(settings)))
        stop("The 'settings' argument does not contain column 'N'.")
    if (!checkmate::test_numeric(settings[, "K"], lower = 1))
        stop("'K' contains a non-positive value.")
    if (!checkmate::test_numeric(settings[, "N"], lower = 1))
        stop("'N' contains a non-positive value.")

    # Assert that fit is an occumbFit object
    assert_occumbFit(fit)

    # Assert that model parameters are not replicate-specific
    if (!length(dim(get_post_samples(fit, "theta"))) < 4)
        stop("'theta' is replicate-specific: the current 'eval_util_L' is not applicable to models with replicate-specific parameters.")
    if (!length(dim(get_post_samples(fit, "phi"))) < 4)
        stop("'phi' is replicate-specific: the current 'eval_util_L' is not applicable to models with replicate-specific parameters.")
}

# @title Monte-Carlo integration to obtain expected utility.
# @param z A species presence-absence array
#   (dim = N_sample * N_species * N_sites).
# @param theta A sequence capture probability array
#   (dim = N_sample * N_species * N_sites).
# @param phi A relative sequence dominance array
#   (dim = N_sample, N_species * N_sites).
# @param K The number of replicates (integer).
# @param N Sequence depth (numeric).
# @param scale Spatial scale to evaluate detection effectiveness.
# @param N_rep Controls the sample size for Monte Carlo simulation.
#   The integral is evaluated using a total of N_sample * N_rep random samples.
# @param core The number of cores to use for parallelization.
# @return The expected utility.
eutil <- function(z, theta, phi, K, N, scale = c("local", "regional"),
                  N_rep = N_rep, cores = cores) {

    M <- dim(z)[1]

    if (match.arg(scale) == "local") {
        fun <- .cutil_local
    } else if (match.arg(scale) == "regional") {
        fun <- .cutil_regional
    }

    if (cores == 1) {
        util_rep <- unlist(
            lapply(X = rep(seq_len(M), each = N_rep),
                   FUN = fun,
                   args = list(z = z, theta = theta, phi = phi, K = K, N = N))
        )
    } else {
        if (.Platform$OS.type == "windows") {
            # On Windows use makePSOCKcluster() and parLapply() for multiple cores
            cl <- parallel::makePSOCKcluster(cores)
            on.exit(parallel::stopCluster(cl))
            util_rep <- unlist(
                parallel::parLapply(cl = cl,
                                    X = rep(seq_len(M), each = N_rep),
                                    FUN = fun,
                                    args = list(z = z,
                                                theta = theta,
                                                phi = phi,
                                                K = K,
                                                N = N))
            )
        } else {
            # On Mac or Linux use mclapply() for multiple cores
            util_rep <- unlist(
                parallel::mclapply(mc.cores = cores,
                                   X = rep(seq_len(M), each = N_rep),
                                   FUN = fun,
                                   args = list(z = z,
                                               theta = theta,
                                               phi = phi,
                                               K = K,
                                               N = N))
            )
        }
    }

    mean(util_rep)
}

# @title Conditional utility function for local species diversity assessments.
# @param z A species presence-absence matrix (dim = N_species * N_sites).
# @param theta A sequence capture probability matrix (dim = N_species * N_sites).
# @param phi A relative sequence dominance matrix (dim = N_species * N_sites).
# @param K The number of replicates (integer).
# @param N Sequence depth (numeric).
# @return The expected number of detected species per site.
cutil_local <- function(z, theta, phi, K, N) {
    pi <- predict_pi(z, theta, phi, K)
    detect_probs <- predict_detect_probs_local(pi, N)
    return(sum(detect_probs) / ncol(z))
}

# @title Conditional utility function for regional species diversity assessments.
# @param z A species presence-absence matrix (dim = N_species * N_sites).
# @param theta A sequence capture probability matrix (dim = N_species * N_sites).
# @param phi A relative sequence dominance matrix (dim = N_species * N_sites).
# @param K The number of replicates (integer).
# @param N Sequence depth (numeric).
# @return The expected total number of species detected over the all sites.
cutil_regional <- function(z, theta, phi, K, N) {
    pi <- predict_pi(z, theta, phi, K)
    detect_probs <- predict_detect_probs_regional(pi, N)
    return(sum(detect_probs))
}

# Auxiliary functions adapting cutil to lapply
.cutil_local <- function(n, args) {
    cutil_local(args$z[n, , ],
                args$theta[n, , ],
                args$phi[n, , ],
                args$K,
                args$N)
}

.cutil_regional <- function(n, args) {
    cutil_regional(args$z[n, , ],
                   args$theta[n, , ],
                   args$phi[n, , ],
                   args$K,
                   args$N)
}

# Calculate pi
predict_pi <- function(z, theta, phi, K) {
    I <- nrow(z)
    J <- ncol(z)

    u <- r <- array(dim = c(I, J, K))
    pi <- array(dim = c(I, J, K))

    count <- 0
    for (n in seq_len(1000)) {
        # Posterior predictive samples of u and r
        for (k in seq_len(K)) {
            u[, , k] <- stats::rbinom(I * J, 1, z * theta)
            r[, , k] <- stats::rgamma(I * J, phi, 1)
        }

        # Derive pi
        for (j in seq_len(J)) {
            for (k in seq_len(K)) {
                ur <- (u * r)[, j, k]
                pi[, j, k] <- ur / sum(ur)
            }
        }

        if (any(is.nan(pi))) count <- count + 1
        else break
    }

    if (count == 1000)
        stop("Failed to generate valid pi values under the given parameter set.")

    pi
}

# Calculate species detection probabilities (local scale)
predict_detect_probs_local <- function(pi, N) {
    I <- dim(pi)[1]
    J <- dim(pi)[2]

    detect_probs <- matrix(nrow = I, ncol = J)
    for (i in seq_len(I)) {
        for (j in seq_len(J)) {
            detect_probs[i, j] <- 1 - prod((1 - pi[i, j, ])^N)
        }
    }

    detect_probs
}

# Calculate species detection probabilities (regional scale)
predict_detect_probs_regional <- function(pi, N) {
    I <- dim(pi)[1]

    detect_probs <- vector(length = I)
    for (i in seq_len(I)) {
        detect_probs[i] <- 1 - prod((1 - pi[i, , ])^N)
    }

    detect_probs
}

