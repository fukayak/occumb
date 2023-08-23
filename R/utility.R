#' @title Expected utility for local species diversity assessments.
#'
#' @description `eval_util_L()` evaluates the expected utility for the local
#'  species diversity assessment using Monte Carlo integration.
#' @details
#' The utility of the local species diversity assessment can be defined, for a
#'  given set of sites, as the expected number of detected species per site
#'  (Fukaya et al. 2022). `eval_util_L()` evaluates this utility for arbitrary
#'  sets of sites that can potentially have different values for site occupancy
#'  status of species, \eqn{z}{z}, sequence capture probabilities of species,
#'  \eqn{\theta}{theta}, and sequence relative dominance of species,
#'  \eqn{\phi}{phi}, for the combination of `K` and `N` values specified in the
#'  `conditions` argument.
#'  Such evaluations can be used to balance `K` and `N` to maximize the utility
#'  under a constant budget (possible combinations of `K` and `N` under
#'  specified budget and cost values are easily obtained using `list_cond_L()`;
#'  see the example below).
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
#' Expectations are taken with respect to posterior (or, possibly prior)
#'  predictive distributions of \eqn{\boldsymbol{r} = \{r_{ijk}\}}{r} and
#'  \eqn{\boldsymbol{u} = \{u_{ijk}\}}{u}, which are evaluated numerically by
#'  Monte Carlo integration. The predictive distributions of
#'  \eqn{\boldsymbol{r}}{r} and \eqn{\boldsymbol{u}}{u} depend on the model
#'  parameters \eqn{z}{z}, \eqn{\theta}{theta}, and \eqn{\phi}{phi} values.
#'  Their posterior (or prior) distribution is specified by supplying an
#'  `occumbFit` object containing their posterior samples via the `fit` argument
#'  or by supplying a matrix or array of posterior (or prior) samples of
#'  parameter values via the `z`, `theta`, and `phi` arguments. Higher
#'  approximation accuracy can be obtained by increasing the value of `N_rep`.
#'
#' The `eval_util_L()` function can be executed by supplying the `fit` argument
#'  without specifying the `z`, `theta`, and `phi` arguments, by supplying the
#'  three `z`, `theta`, and `phi` arguments without the `fit` argument, or by
#'  supplying the `fit` argument and any or all of the `z`, `theta`, and `phi`
#'  arguments.
#'
#' The expected utility is evaluated assuming homogeneity of replicates, in the
#'  sense that \eqn{\theta}{theta} and \eqn{\phi}{phi}, the model parameters
#'  associated with the species detection process, are constant across
#'  replicates within a site. For this reason, `eval_util_L()` does not accept
#'  replicate-specific \eqn{\theta}{theta} and \eqn{\phi}{phi}. If the
#'  `occumbFit` object supplied in the `fit` argument has a replicate-specific
#'  parameter, parameter samples to be used in the utility evaluation must be
#'  provided explicitly via `theta` or `phi` argument.
#'
#' Monte Carlo integration is executed in parallel on multiple CPU cores where
#'  the `cores` argument controls the degree of parallelization.
#' @param settings A data frame that specifies a set of conditions under which
#'  the utility is evaluated. It must include a column named `K` and `N`, which
#'  specifies the number of replicates per site and the sequencing depth per
#'  replicate, respectively.
#'  `K` and `N` must be numeric vectors greater than 0. When `K` contains a
#'  decimal, the decimal part is discarded and treated as an integer.
#'  Additional columns are ignored but may be included.
#' @param fit An `occumbFit` object containing a posterior sample of the
#'  relevant parameters.
#' @param z Sample values of site occupancy status of species stored in an array
#'  with sample \eqn{\times}{*} species \eqn{\times}{*} site dimension.
#' @param theta Sample values of sequence capture probabilities of species
#'  stored in a matrix with sample \eqn{\times}{*} species dimension or an array
#'  with sample \eqn{\times}{*} species \eqn{\times}{*} site dimension.
#' @param phi Sample values of sequence relative dominance of species stored in
#'  a matrix with sample \eqn{\times}{*} species dimension or an array with
#' sample \eqn{\times}{*} species \eqn{\times}{*} site dimension.
#' @param N_rep Controls the sample size for Monte Carlo integration.
#'  The integral is evaluated using a total of `N_sample * N_rep` random samples,
#'  where `N_sample` is the maximum size of the MCMC sample in the `fit`
#'  argument and the parameter sample in the `z`, `theta`, and `phi` arguments.
#' @param cores The number of cores to use for parallelization.
#' @return A data frame with a column named `Utility` in which the estimates of
#'  expected utility are stored. This is obtained by adding the `Utility` column
#'  to the data frame provided in the `settings` argument.
#' @section References:
#'      K. Fukaya, N. I. Kondo, S. S. Matsuzaki and T. Kadoya (2022)
#'      Multispecies site occupancy modelling and study design for spatially
#'      replicated environmental DNA metabarcoding. \emph{Methods in Ecology
#'      and Evolution} \strong{13}:183--193.
#'      \doi{10.1111/2041-210X.13732}
#' @examples
#' \dontrun{
#' # Generate the smallest random dataset (2 species * 2 sites * 2 reps)
#' I <- 2 # Number of species
#' J <- 2 # Number of sites
#' K <- 2 # Number of replicates
#' data <- occumbData(
#'     y = array(sample.int(I * J * K), dim = c(I, J, K)),
#'     spec_cov = list(cov1 = rnorm(I)),
#'     site_cov = list(cov2 = rnorm(J),
#'                     cov3 = factor(1:J)),
#'     repl_cov = list(cov4 = matrix(rnorm(J * K), J, K)))
#'
#' # Fitting a null model
#' fit <- occumb(data = data)
#'
#' # Estimate expected utility
#' util1 <- eval_util_L(expand.grid(K = 1:3, N = c(1E3, 1E4, 1E5)),
#'                      fit) # Arbitrary K and N values
#' util2 <- eval_util_L(list_cond_L(budget = 1E5,
#'                                  lambda1 = 0.01,
#'                                  lambda2 = 5000,
#'                                  fit),
#'                      fit) # K and N values under specified budget and cost
#' util3 <- eval_util_L(list_cond_L(budget = 1E5,
#'                                  lambda1 = 0.01,
#'                                  lambda2 = 5000,
#'                                  fit,
#'                                  K = 1:5),
#'                      fit) # K values restricted
#' }
#' @export
eval_util_L <- function(settings,
                        fit = NULL,
                        z = NULL,
                        theta = NULL,
                        phi = NULL,
                        N_rep = 1,
                        cores = 1L) {

    # Validate arguments
    qc_eval_util_L(settings, fit, z, theta, phi)

    # Set parameter values
    if (is.null(z))
        z <- get_post_samples(fit, "z")

    if (is.null(theta))
        theta <- get_post_samples(fit, "theta")

    if (is.null(phi))
        phi <- get_post_samples(fit, "phi")

    # Resampling to match sample size
    n_samples <- c(dim(z)[1], dim(theta)[1], dim(phi)[1])

    if (n_samples[1] < max(n_samples))
        z <- z[sample.int(n_samples[1], max(n_samples), replace = TRUE), , ]

    if (n_samples[2] < max(n_samples)) {
        if (length(dim(theta)) == 2)
            theta <- theta[sample.int(n_samples[2], max(n_samples), replace = TRUE), ]
        if (length(dim(theta)) == 3)
            theta <- theta[sample.int(n_samples[2], max(n_samples), replace = TRUE), , ]
    }

    if (n_samples[3] < max(n_samples)) {
        if (length(dim(phi)) == 2)
            phi <- phi[sample.int(n_samples[3], max(n_samples), replace = TRUE), ]
        if (length(dim(phi)) == 3)
            phi <- phi[sample.int(n_samples[3], max(n_samples), replace = TRUE), , ]
    }

    # Adapt theta/phi when they are species-specific
    if (length(dim(theta)) == 2)
        theta <- outer(theta, rep(1, dim(z)[3]))
    if (length(dim(phi)) == 2)
        phi <- outer(phi, rep(1, dim(z)[3]))

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

#' @title Expected utility for regional species diversity assessments.
#'
#' @description `eval_util_R()` evaluates the expected utility for the regional
#'  species diversity assessment using Monte Carlo integration.
#' @details
#' The utility of the regional species diversity assessment can be defined as
#'  the number of species in the region of interest expected to be detected
#'  (Fukaya et al. 2022). `eval_util_R()` evaluates this utility for the region
#'  modeled in the `occumbFit` object for the combination of `J`, `K`, and `N`
#'  values specified in the `conditions` argument.
#'  Such evaluations can be used to balance `J`, `K`, and `N` to maximize the
#'  utility under a constant budget (possible combinations of `J`, `K`, and `N`
#'  under specified budget and cost values are easily obtained using
#'  `list_cond_R()`; see the example below).
#'  It is also possible to examine how the utility varies with different `J`,
#'  `K`, and `N` values without setting a budget level, which may be useful to know
#'  the satisfactory level of `J`, `K`, and `N` from a purely technical point of
#'  view.
#' 
#' The expected utility is defined as the expected value of the conditional
#'  utility of the form:
#'  \deqn{U(J, K, N \mid \boldsymbol{r}, \boldsymbol{u}) = \sum_{i = 1}^{I}\left\{1 - \prod_{j = 1}^{J}\prod_{k = 1}^{K}\left(1 - \frac{u_{ijk}r_{ijk}}{\sum_{m = 1}^{I}u_{mjk}r_{mjk}} \right)^N \right\}}{U(J, K, N | r, u) = sum_{i}((1 - \prod_{j}\prod_{k}(1 - (u[i, j, k] * r[i, j, k])/sum(u[, j, k] * r[, j, k])))^N)}
#'  where \eqn{u_{ijk}}{u[i, j, k]} is a latent indicator variable representing
#'  the inclusion of the sequence of species \eqn{i}{i} in replicate \eqn{k}{k}
#'  at site \eqn{j}{j} and \eqn{r_{ijk}}{r[i, j, k]} is a latent variable that
#'  is proportional to the relative frequency of the sequence of species
#'  \eqn{i}{i}, conditional on its presence in replicate \eqn{k}{k} at site
#'  \eqn{j}{j} (Fukaya et al. 2022).
#' Expectations are taken with respect to the posterior predictive distributions
#'  of \eqn{\boldsymbol{r} = \{r_{ijk}\}}{r} and
#'  \eqn{\boldsymbol{u} = \{u_{ijk}\}}{u}, which are evaluated numerically by
#'  Monte Carlo integration using MCMC samples in the `occumbFit` object. Higher
#'  approximation accuracy can be obtained by increasing the value of `N_rep`.
#'
#' The expected utility is evaluated assuming that all replicates are
#'  homogeneous in the sense that the model parameters are constant across
#'  replicates. For this reason, in the current version, `eval_util_R()`
#'  cannot be applied if the supplied `occumbFit` object contains a model with
#'  replicate covariates.
#'
#' If the relevant parameters are modeled as a function of site covariates, the
#'  expected utility is evaluated to account for the site heterogeneity of the
#'  parameters. Specifically, to incorporate site heterogeneity, the parameter
#'  values for each `J` site are determined by selecting site-specific
#'  parameters estimated in the model via random sampling with replacement.
#'  Thus, the expected utility is evaluated by assuming the set of parameter
#'  values for the sites modeled in the `fit` as a statistical population of
#'  site-specific parameters.
#'
#' Monte Carlo integration is executed in parallel on multiple CPU cores where
#'  the `cores` argument controls the degree of parallelization.
#' @param settings A data frame that specifies a set of conditions under which
#'  the utility is evaluated. It must include a column named `J`, `K`, and `N`,
#'  which specifies the number of sites, the number of replicates per site, and
#'  the sequencing depth per replicate, respectively.
#'  `J`, `K`, and `N` must be numeric vectors greater than 0. When `J` and `K`
#'  contains a decimal, the decimal part is discarded and treated as an integer.
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
#'      \doi{10.1111/2041-210X.13732}
#' @examples
#' \dontrun{
#' # Generate the smallest random dataset (2 species * 2 sites * 2 reps)
#' I <- 2 # Number of species
#' J <- 2 # Number of sites
#' K <- 2 # Number of replicates
#' data <- occumbData(
#'     y = array(sample.int(I * J * K), dim = c(I, J, K)),
#'     spec_cov = list(cov1 = rnorm(I)),
#'     site_cov = list(cov2 = rnorm(J),
#'                     cov3 = factor(1:J)),
#'     repl_cov = list(cov4 = matrix(rnorm(J * K), J, K)))
#'
#' # Fitting a null model
#' fit <- occumb(data = data)
#'
#' # Estimate expected utility
#' util1 <- eval_util_R(expand.grid(J = 1:3, K = 1:3, N = c(1E3, 1E4, 1E5)),
#'                      fit) # Arbitrary J, K, and N values
#' util2 <- eval_util_R(list_cond_R(budget = 50000,
#'                                  lambda1 = 0.01,
#'                                  lambda2 = 5000,
#'                                  lambda3 = 5000),
#'                      fit) # J, K, and N values under specified budget and cost
#' util3 <- eval_util_R(list_cond_R(budget = 50000,
#'                                  lambda1 = 0.01,
#'                                  lambda2 = 5000,
#'                                  lambda3 = 5000,
#'                                  K = 1:5),
#'                      fit) # K values restricted
#' util4 <- eval_util_R(list_cond_R(budget = 50000,
#'                                  lambda1 = 0.01,
#'                                  lambda2 = 5000,
#'                                  lambda3 = 5000,
#'                                  J = 1:3, K = 1:5),
#'                      fit) # J and K values restricted
#' }
#' @export
eval_util_R <- function(settings,
                        fit,
                        N_rep = 1,
                        cores = 1L) {

    # Validate arguments
    qc_eval_util_R(settings, fit)

    # Extract posterior samples
    psi   <- get_post_samples(fit, "psi")
    theta <- get_post_samples(fit, "theta")
    phi   <- get_post_samples(fit, "phi")

    # Determine site dimension
    has_site_dim <- c(length(dim(psi))   == 3,
                      length(dim(theta)) == 3,
                      length(dim(phi))   == 3)
    N_site <-
        if (has_site_dim[1]) {
            dim(psi)[3]
        } else if (has_site_dim[2]) {
            dim(theta)[3]
        } else if (has_site_dim[3]) {
            dim(phi)[3]
        }

    # Make N_rep copies of psi, theta, and phi
    if (N_rep > 1) {
        psi   <- copy_psi_theta_phi(psi, has_site_dim[1], N_rep)
        theta <- copy_psi_theta_phi(theta, has_site_dim[2], N_rep)
        phi   <- copy_psi_theta_phi(phi, has_site_dim[3], N_rep)
    }

    # Calculate expected utility
    result <- rep(NA, nrow(settings))
    for (i in seq_len(nrow(settings))) {
        # Random sampling of sites
        J <- settings[i, "J"]
        site_use <- matrix(nrow = dim(psi)[1], ncol = J)
        if (any(has_site_dim)) {
            for (n in seq_len(nrow(site_use)))
                site_use[n, ] <- sample.int(N_site, J, replace = TRUE)
        }

        # Prepare z, theta, and phi (dim = N_sample * N_species * N_site)
        z_i     <- get_z_i(psi, has_site_dim[1], site_use)
        theta_i <- get_theta_phi_i(theta, has_site_dim[2], site_use)
        phi_i   <- get_theta_phi_i(phi, has_site_dim[3], site_use)

        result[i] <- eutil(z = z_i, theta = theta_i, phi = phi_i,
                           K = settings[i, "K"], N = settings[i, "N"],
                           scale = "regional",
                           N_rep = 1, cores = cores)
    }

    # Output
    out <- cbind(settings, result)
    colnames(out)[ncol(out)] <- "Utility"
    out
}

#' @title Conditions for local assessment under certain budget and cost values.
#' @description `list_cond_L()` constructs a list of possible local species
#'  diversity assessment conditions under specified budget and cost values.
#' @details
#'   This function can generate a data frame object to be given to the
#'  `settings` argument of `eval_util_L()`; see Examples of `eval_util_L()`.
#'  By default, it outputs a list of all feasible combinations of values for the
#'  number of replicates per site `K` and the sequencing depth per replicate
#'  `N`, based on the given budget and cost values and the number of sites
#'  (identified by reference to the `fit` object). The resulting `N` can be
#'  non-integer because it is calculated simply by assuming that we can obtain
#'  its maximum value. If you want to obtain a list for only a subset of the
#'  possible `K` values under a given budget and cost values, use the `K`
#'  argument to provide a vector of the desired `K` values.
#' @param budget A numeric specifying the amount of budget. Use the currency
#'  unit consistent with `lambda1` and `lambda2`.
#' @param lambda1 A numeric specifying the cost per sequence read for
#'  high-throughput sequencing. Use the currency unit consistent with `budget`
#'  and `lambda2`.
#' @param lambda2 A numeric specifying the cost per replicate for library
#'  preparation. Use the currency unit consistent with `budget` and `lambda1`.
#' @param fit An `occumbFit` object.
#' @param K An optional vector for manually specifying the number of replicates.
#' @return A data frame containing columns named `budget`, `lambda1`, `lambda2`,
#'  `K`, and `N`.
#' @export
list_cond_L <- function(budget, lambda1, lambda2, fit, K = NULL) {

    ## Validate arguments
    # Assert that budget and cost values are positive.
    if (!checkmate::test_numeric(budget, lower = 0))
        stop("Negative 'budget' value.\n")
    if (!checkmate::test_numeric(lambda1, lower = 0))
        stop("Negative 'lambda1' value.\n")
    if (!checkmate::test_numeric(lambda2, lower = 0))
        stop("Negative 'lambda2' value.\n")

    # Assert that fit is an occumbFit object
    assert_occumbFit(fit)

    # Determine the number of sites
    J <- dim(get_post_samples(fit, "z"))[3]

    # Determine max_K under given budget, cost, and the number of sites
    max_K <- floor(budget / (lambda2 * J))

    # Assert that given settings ensure at least one replicate per site
    if (!checkmate::test_numeric(max_K, lower = 1))
        stop("Impossible to have > 0 replicates per site under the given budget, cost, and the number of sites.\n")

    if (is.null(K)) {
        # Determine K values
        K <- seq_len(max_K)[budget - lambda2 * J * seq_len(max_K) > 0]
    } else {
        # Assert that K >= 1
        if (!checkmate::test_numeric(K, lower = 1))
            stop("'K' contains values less than one.\n")

        # Assert that the all given K are feasible
        if (any(!budget - lambda2 * J * K > 0))
            stop(paste("A value of 'K' greater than",
                       max_K,
                       "is not feasible under the given budget, cost, and the number of sites.\n"))
    }

    ## Output a table of conditions
    out <- cbind(rep(budget, length(K)),
                 rep(lambda1, length(K)),
                 rep(lambda2, length(K)),
                 K,
                 (budget - lambda2 * J * K) / (lambda1 * J * K))
    colnames(out) <- c("budget", "lambda1", "lambda2", "K", "N")
    data.frame(out)
}

#' @title Conditions for regional assessment under certain budget and cost values.
#' @description `list_cond_R()` constructs a list of possible regional species
#'  diversity assessment conditions under specified budget and cost values.
#' @details
#'   This function can generate a data frame object to be given to the
#'  `settings` argument of `eval_util_R()`; see Examples of `eval_util_R()`.
#'  By default, it outputs a list of all feasible combinations of values for the
#'  number of sites `J`, number of replicates per site `K`, and the sequencing
#'  depth per replicate `N`, based on the given budget and cost values.
#'  The resulting `N` can be non-integer because it is calculated simply by
#'  assuming that we can obtain its maximum value.
#'  If you want to obtain a list for only a subset of the possible values of `J`
#'  and `K` under a given budget and cost values, use the `J` and/or `K`
#'  arguments (in fact, it is recommended that a relatively small number of `K`
#'  values be specified using the `K` argument because the list of all
#'  conditions achievable under moderate budget and cost values can be huge, and
#'  it is rarely practical to have a vast number of replicates per site). If a
#'  given combination of `J` and `K` values is not feasible under the specified
#'  budget and cost values, that combination will be ignored and excluded from
#'  the output.
#' @param budget A numeric specifying the amount of budget. Use the currency
#'  unit consistent with `lambda1`, `lambda2`, and `lambda3`.
#' @param lambda1 A numeric specifying the cost per sequence read for
#'  high-throughput sequencing. Use the currency unit consistent with `budget`,
#'  `lambda2`, and `lambda3`.
#' @param lambda2 A numeric specifying the cost per replicate for library
#'  preparation. Use the currency unit consistent with `budget`, `lambda1`, and
#'  `lambda3`.
#' @param lambda3 A numeric specifying the visiting cost per site. Use the
#'  currency unit consistent with `budget`, `lambda1`, and `lambda2`.
#' @param J An optional vector for manually specifying the number of sites
#' @param K An optional vector for manually specifying the number of replicates.
#'  For computational convenience, `K` values must be in ascending order.
#' @return A data frame containing columns named `budget`, `lambda1`, `lambda2`,
#'  , `lambda3`, `J`, `K`, and `N`.
#' @export
list_cond_R <- function(budget, lambda1, lambda2, lambda3, J = NULL, K = NULL) {

    ## Validate arguments
    # Assert that budget and cost values are positive.
    if (!checkmate::test_numeric(budget, lower = 0))
        stop("Negative 'budget' value.\n")
    if (!checkmate::test_numeric(lambda1, lower = 0))
        stop("Negative 'lambda1' value.\n")
    if (!checkmate::test_numeric(lambda2, lower = 0))
        stop("Negative 'lambda2' value.\n")
    if (!checkmate::test_numeric(lambda3, lower = 0))
        stop("Negative 'lambda3' value.\n")

    # Assert that J, K >= 1
    if (!is.null(J) & !checkmate::test_numeric(J, lower = 1))
        stop("'J' contains values less than one.\n")
    if (!is.null(K) & !checkmate::test_numeric(K, lower = 1))
        stop("'K' contains values less than one.\n")

    # Assert that K is in ascending order
    if (!is.null(K) & !identical(K, sort(K)))
        stop("'K' must be in ascending order.\n")

    # Determine the combination of J and K to be used
    if (is.null(J)) {
        max_J <- find_maxJ(budget, lambda2, lambda3)
        if (!max_J)
            stop("No valid combination of 'J' and 'K' under the given budget and cost.\n")
        J <- seq_len(max_J)
    }
    if (is.null(K)) {
        max_K <- find_maxK(budget, lambda2, lambda3)
        if (!max_K)
            stop("No valid combination of 'J' and 'K' under the given budget and cost.\n")
        K <- seq_len(max_K)
    }

    J_valid <- K_valid <- vector()
    for (k in K) {
        J_valid_k <- J[budget - lambda2 * J * k - lambda3 * J > 0]
        if (length(J_valid_k) > 0) {
            J_valid <- c(J_valid, J_valid_k)
            K_valid <- c(K_valid, rep(k, length(J_valid_k)))
        } else {
            break
        }
    }

    # Assert that given settings ensure at least one valid combination of J and K
    if (!length(J_valid) > 0)
        stop("No valid combination of 'J' and 'K' under the given budget and cost.\n")

    ## Output a table of conditions
    out <- cbind(rep(budget, length(J_valid)),
                 rep(lambda1, length(J_valid)),
                 rep(lambda2, length(J_valid)),
                 rep(lambda3, length(J_valid)),
                 J_valid,
                 K_valid,
                 (budget - lambda2 * J_valid * K_valid - lambda3 * J_valid) / (lambda1 * J_valid * K_valid))
    colnames(out) <- c("budget", "lambda1", "lambda2", "lambda3", "J", "K", "N")
    data.frame(out)
}

qc_eval_util_L <- function(settings, fit, z, theta, phi) {
    # Assert that settings is a data frame and contains the required columns
    checkmate::assert_data_frame(settings)
    if (!checkmate::testSubset("K", names(settings)))
        stop("The 'settings' argument does not contain column 'K'.\n")
    if (!checkmate::testSubset("N", names(settings)))
        stop("The 'settings' argument does not contain column 'N'.\n")
    if (!checkmate::test_numeric(settings[, "K"], lower = 1))
        stop("'K' contains values less than one.\n")
    if (!checkmate::test_numeric(settings[, "N"], lower = 1))
        stop("'N' contains values less than one.\n")

    # Assert that either fit or (z, theta, phi) is provided
    if (is.null(fit) & !(!is.null(z) & !is.null(theta) & !is.null(phi)))
        stop("Parameter values are not fully specified: use fit argument or otherwise use all of z, theta, phi arguments.\n")

    if (!is.null(fit)) {
        # Assert that fit is an occumbFit object
        assert_occumbFit(fit)

        # Stop when modeled parameters are replicate-specific
        if (length(dim(get_post_samples(fit, "theta"))) == 4)
            stop("'fit' contains replicate-specific theta: specify appropriate theta values via the 'theta' argument to run.\n")
        if (length(dim(get_post_samples(fit, "phi"))) == 4)
            stop("'fit' contains replicate-specific phi: specify appropriate phi values via the 'phi' argument to run.\n")
    }

    # Assert that z, theta, and phi have an appropriate dimension and elements
    checkmate::assert_array(z, d = 3, null.ok = TRUE)
    checkmate::assert_array(theta, min.d = 2, max.d = 3, null.ok = TRUE)
    checkmate::assert_array(phi, min.d = 2, max.d = 3, null.ok = TRUE)
    checkmate::assert_integerish(z, lower = 0, upper = 1,
                                 any.missing = FALSE, null.ok = TRUE)
    checkmate::assert_numeric(theta, lower = 0, upper = 1,
                              any.missing = FALSE, null.ok = TRUE)
    checkmate::assert_numeric(phi, lower = 0, any.missing = FALSE, null.ok = TRUE)

    # Assert that dimensions of z, theta, and phi agree with fit
    if (!is.null(fit)) {
        I <- dim(fit@data@y)[1]
        J <- dim(fit@data@y)[2]

        if (!is.null(z)) {
            if (dim(z)[2] != I)
                stop(paste0("Mismatch in species dimension: dim(z)[2] must be ", I, ".\n"))
            if (dim(z)[3] != J)
                stop(paste0("Mismatch in site dimension: dim(z)[3] must be ", J, ".\n"))
        }
        if (!is.null(theta)) {
            if (dim(theta)[2] != I)
                stop(paste0("Mismatch in species dimension: dim(theta)[2] must be ", I, ".\n"))
            if (length(dim(theta)) == 3 & dim(theta)[3] != J)
                stop(paste0("Mismatch in site dimension: dim(theta)[3] must be ", J, ".\n"))
        }
        if (!is.null(phi)) {
            if (dim(phi)[2] != I)
                stop(paste0("Mismatch in species dimension: dim(phi)[2] must be ", I, ".\n"))
            if (length(dim(phi)) == 3 & dim(phi)[3] != J)
                stop(paste0("Mismatch in site dimension: dim(phi)[3] must be ", J, ".\n"))
        }
    } else {
        has_site_dim <- c(TRUE, length(dim(theta)) == 3, length(dim(phi)) == 3)
        vI <- c(dim(z)[2], dim(theta)[2], dim(phi)[2])
        vJ <- c(dim(z)[3],
                ifelse(has_site_dim[2], dim(theta)[3], NA),
                ifelse(has_site_dim[3], dim(phi)[3], NA))
        vJ <- vJ[has_site_dim]
        terms <- c("dim(z)[3]", "dim(theta)[3]", "dim(phi)[3]")[has_site_dim]

        if (length(unique(vI)) != 1)
            stop("Mismatch in species dimension: dim(z)[2], dim(theta)[2], and dim(phi)[2] must be equal.\n")
        if (length(unique(vJ)) != 1)
            stop(paste0("Mismatch in site dimension: ",
                        knitr::combine_words(terms), " must be equal.\n"))
    }

    invisible(NULL)
}

qc_eval_util_R <- function(settings, fit) {
    # Assert that settings is a data frame and contains the required columns
    checkmate::assert_data_frame(settings)
    if (!checkmate::testSubset("J", names(settings)))
        stop("The 'settings' argument does not contain column 'J'.\n")
    if (!checkmate::testSubset("K", names(settings)))
        stop("The 'settings' argument does not contain column 'K'.\n")
    if (!checkmate::testSubset("N", names(settings)))
        stop("The 'settings' argument does not contain column 'N'.\n")
    if (!checkmate::test_numeric(settings[, "J"], lower = 1))
        stop("'J' contains values less than one.\n")
    if (!checkmate::test_numeric(settings[, "K"], lower = 1))
        stop("'K' contains values less than one.\n")
    if (!checkmate::test_numeric(settings[, "N"], lower = 1))
        stop("'N' contains values less than one.\n")

    # Assert that fit is an occumbFit object
    assert_occumbFit(fit)

    # Assert that model parameters are not site- or replicate-specific
    if (length(dim(get_post_samples(fit, "theta"))) == 4)
        stop("'theta' is replicate-specific: the current 'eval_util_R' is not applicable to models with replicate-specific parameters.\n")
    if (length(dim(get_post_samples(fit, "phi"))) == 4)
        stop("'phi' is replicate-specific: the current 'eval_util_R' is not applicable to models with replicate-specific parameters.\n")
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
            parallel::clusterEvalQ(cl, library(occumb))
            on.exit(parallel::stopCluster(cl))
            util_rep <- unlist(
                parallel::parLapply(cl = cl,
                                    X = rep(seq_len(M), each = N_rep),
                                    fun = fun,
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
    if (dim(args$z)[2] == 1) {          # When I = 1
        cutil_local(matrix(args$z[n, , ], nrow = 1),
                    matrix(args$theta[n, , ], nrow = 1),
                    matrix(args$phi[n, , ], nrow = 1),
                    args$K,
                    args$N)
    } else if (dim(args$z)[3] == 1) {   # When J = 1
        cutil_local(matrix(args$z[n, , ], ncol = 1),
                    matrix(args$theta[n, , ], ncol = 1),
                    matrix(args$phi[n, , ], ncol = 1),
                    args$K,
                    args$N)
    } else {
        cutil_local(args$z[n, , ],
                    args$theta[n, , ],
                    args$phi[n, , ],
                    args$K,
                    args$N)
    }
}

.cutil_regional <- function(n, args) {
    if (dim(args$z)[2] == 1) {          # When I = 1
        cutil_regional(matrix(args$z[n, , ], nrow = 1),
                       matrix(args$theta[n, , ], nrow = 1),
                       matrix(args$phi[n, , ], nrow = 1),
                       args$K,
                       args$N)
    } else if (dim(args$z)[3] == 1) {   # When J = 1
        cutil_regional(matrix(args$z[n, , ], ncol = 1),
                       matrix(args$theta[n, , ], ncol = 1),
                       matrix(args$phi[n, , ], ncol = 1),
                       args$K,
                       args$N)
    } else {
        cutil_regional(args$z[n, , ],
                       args$theta[n, , ],
                       args$phi[n, , ],
                       args$K,
                       args$N)
    }
}

# Calculate pi
predict_pi <- function(z, theta, phi, K) {
    I <- nrow(z)
    J <- ncol(z)

    u <- r <- array(dim = c(I, J, K))
    pi <- array(dim = c(I, J, K))

    # Posterior predictive samples of u and r
    for (k in seq_len(K)) {
        u[, , k] <- sample_u(z * theta)
        r[, , k] <- stats::rgamma(I * J, phi, 1)
    }

    # Derive pi
    for (j in seq_len(J)) {
        for (k in seq_len(K)) {
            ur <- (u * r)[, j, k]
            pi[, j, k] <- ur / sum(ur)
        }
    }

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

# Find the maximum K value under the specified budget for regional assessment
find_maxK <- function(budget, lambda2, lambda3, ulim = 1E4) {
    J <- 1
    for (n in 2:log10(ulim)) {
        K <- seq_len(10^n)

        if (any(budget - lambda2 * J * K - lambda3 * J > 0)) {
            maxK <- max(K[budget - lambda2 * J * K - lambda3 * J > 0])
        } else {
            maxK <- 0
        }

        if (maxK < length(K))
            break

        if (n == log10(ulim))
            stop("Maximum `K` value seems too large under the specified budget and cost values: consider using the `K` argument to specify a smaller set of `K` values of interest.\n")
    }
    maxK
}

# Find the maximum J value under the specified budget for regional assessment
find_maxJ <- function(budget, lambda2, lambda3, ulim = 1E6) {
    K <- 1
    for (n in 4:log10(ulim)) {
        J <- seq_len(10^n)

        if (any(budget - lambda2 * J * K - lambda3 * J > 0)) {
            maxJ <- max(J[budget - lambda2 * J * K - lambda3 * J > 0])
        } else {
            maxJ <- 0
        }

        if (maxJ < length(J))
            break

        if (n == log10(ulim))
            stop("Maximum `J` value seems too large under the specified budget and cost values: consider using the `J` argument to specify a smaller set of `J` values of interest.\n")
    }
    maxJ
}

# Sample z except for all zero cases
sample_z <- function(psi) {
    z <- stats::rbinom(length(psi), 1, psi)

    # Resample when all species have z = 0
    if (sum(z)) {
        return(z)
    } else {
        while(sum(z) == 0)
            z <- stats::rbinom(length(psi), 1, psi)
        return(z)
    }
}

# Sample u except for all zero cases
sample_u <- function(z_theta) {
    u <- array(stats::rbinom(nrow(z_theta) * ncol(z_theta), 1, z_theta),
               dim = dim(z_theta))

    # Find cases where all species have u = 0 to resample
    allzero <- which(colSums(u) == 0)
    if (length(allzero)) {
        for (n in seq_along(allzero)) {
            while(sum(u[, allzero[n]]) == 0)
                u[, allzero[n]] <-
                    stats::rbinom(nrow(z_theta), 1, z_theta[, allzero[n]])
        }
        return(u)
    } else {
        return(u)
    }
}

get_z_i <- function(psi, has_site_dim, site_use) {
    out <- array(NA, dim = c(dim(psi)[1], dim(psi)[2], ncol(site_use)))

    if (has_site_dim) {
        for (j in seq_len(dim(out)[3])) {
            for (n in seq_len(dim(out)[1]))
                out[n, , j] <- sample_z(psi[n, , site_use[n, j]])
        }
    } else {
        for (j in seq_len(dim(out)[3])) {
            for (n in seq_len(dim(out)[1]))
                out[n, , j] <- sample_z(psi[n, ])
        }
    }

    out
}

get_theta_phi_i <- function(param, has_site_dim, site_use) {
    if (has_site_dim) {
        out <- array(NA, dim = c(dim(param)[1], dim(param)[2], ncol(site_use)))

        for (j in seq_len(dim(out)[3])) {
            for (n in seq_len(dim(out)[1])) {
                out[n, , j] <- param[n, , site_use[n, j]]
            }
        }
    } else {
        out <- outer(param, rep(1, ncol(site_use)))
    }

    out
}

copy_psi_theta_phi <- function(param, has_site_dim, N_rep) {
    if (has_site_dim) {
        out <- param[rep(seq_len(dim(param)[1]), N_rep), , ]
    } else {
        out <- param[rep(seq_len(dim(param)[1]), N_rep), ]
    }
    out
}

