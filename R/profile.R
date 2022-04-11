#' @title Conditional utility function for local species diversity assessments.
#' @param z A species presence-absence matrix (dim = N_species * N_sites).
#' @param theta A sequence capture probability matrix (dim = N_species * N_sites).
#' @param phi A relative sequence dominance matrix (dim = N_species * N_sites).
#' @param K The number of replicates (integer).
#' @param N Sequence depth (numeric).
#' @return The expected number of detected species per site.
cutil_local <- function(z, theta, phi, K, N) {
    pi <- predict_pi(z, theta, phi, K)
    detect_probs <- predict_detect_probs_local(pi, N)
    return(sum(detect_probs) / J)
}

#' @title Conditional utility function for regional species diversity assessments.
#' @param z A species presence-absence matrix (dim = N_species * N_sites).
#' @param theta A sequence capture probability matrix (dim = N_species * N_sites).
#' @param phi A relative sequence dominance matrix (dim = N_species * N_sites).
#' @param K The number of replicates (integer).
#' @param N Sequence depth (numeric).
#' @return The expected total number of species detected over the all sites.
cutil_regional <- function(z, theta, phi, K, N) {
    pi <- predict_pi(z, theta, phi, K)
    detect_probs <- predict_detect_probs_regional(pi, N)
    return(sum(detect_probs))
}

# Calculate pi
predict_pi <- function(z, theta, phi, K) {
    I <- nrow(z)
    J <- ncol(z)

    # Posterior predictive samples of u and r
    u <- r <- array(dim = c(I, J, K))
    for (k in seq_len(K)) {
        u[, , k] <- rbinom(I * J, 1, z * theta)
        r[, , k] <- rgamma(I * J, phi, 1)
    }

    # Derive pi
    pi <- array(dim = c(I, J, K))
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

