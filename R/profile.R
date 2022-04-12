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
# @param rep Controls the sample size for Monte Carlo simulation.
#   The integral is evaluated using a total of N_sample * rep random samples.
# @return The expected utility.
eutil <- function(z, theta, phi, K, N, scale = c("local", "regional"), rep = 1) {

    M <- dim(z)[1]
    util_rep <- vector(length = M * rep)

    if (match.arg(scale) == "local") {
        for (n in rep(seq_len(M), each = rep))
            util_rep[n] <-
                cutil_local(z[n, , ], theta[n, , ], phi[n, , ], K, N)
    }

    if (match.arg(scale) == "regional") {
        for (n in rep(seq_len(M), each = rep))
            util_rep[n] <-
                cutil_regional(z[n, , ], theta[n, , ], phi[n, , ], K, N)
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

