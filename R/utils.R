# Not-in
`%!in%` <- Negate(`%in%`)

# Validation
qc_occumbFit <- function(fit) {
    # Check object classes
    if (!inherits(fit, "occumbFit")) {
        stop("An occumbFit class object is expected for 'fit'")
    } else {
        invisible(NULL)
    }
}

# Multinomial log-likelihood
llmulti <- function(y, N, pi) {
    stats::dmultinom(y, N, pi, log = TRUE)
}

