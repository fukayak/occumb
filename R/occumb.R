#' Model-fitting function.
#' 
#' @param data occumbData class object.
#' @export
occumb <- function(data) {
    # QC

    # Set constants
    const <- set_const(data)

    # Parse the model expression

    # Determine the order of the species effect

    # Set initial values

    # Set parameters monitored

    # Run MCMC in JAGS

    # Embed result in a model object class
}

# Extract constants for the model
set_const <- function(data) {
    y <- data@y
    I <- dim(y)[1]              # Number of species
    J <- dim(y)[2]              # Number of sites
    K <- dim(y)[3]              # Number of replicates
    N <- apply(y, c(2, 3), sum) # Sequence depth

    for (j in 1:J) {
        for (k in 1:K) {
            # Set y to NA and N to 1 when sequence depth is zero.
            if (N[j, k] == 0) {
                y[, j, k] <- NA
                N[j, k]   <- 1
            }
        }
    }

    out <- list(y = y, I = I, J = J, K = K, N = N)
    out
}

