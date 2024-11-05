# Not-in
`%!in%` <- Negate(`%in%`)

# Validation
assert_occumbFit <- function(fit) {
  # Check object classes
  if (!inherits(fit, "occumbFit")) {
    stop("An occumbFit class object is expected for 'fit'\n")
  } else {
    invisible(NULL)
  }
}

# Multinomial log-likelihood
llmulti <- function(y, N, pi) {
  stats::dmultinom(y, N, pi, log = TRUE)
}

# Running a function with a specified seed (used for unit testing)
with_seed <- function(seed, code) {
  code <- substitute(code)
  orig.seed <- .Random.seed
  on.exit(.Random.seed <<- orig.seed)
  set.seed(seed)
  eval.parent(code)
}
