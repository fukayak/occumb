# ------------------------------------------------------------------------------
# NOTE ON ORIGIN / LICENSE (IMPORTANT)
#
# This file is a modified version of code distributed with the book companion
# materials for:
#   "Bayesian population analysis using WinBUGS: A Hierarchical Perspective"
#   (Marc Kéry and Michael Schaub)
#
# Source (book companion website, Swiss Ornithological Institute):
#   https://www.vogelwarte.ch/en/research/population-biology/book-bpa/
#
# In particular, this implementation is adapted from `nimbleSummary.R` included
# in the R package "BPAbook" (e.g., BPAbook_0.0.1.tar.gz) available from the
# above site. The companion site notes that the `nimbleSummary` function was
# written by Ken Kellner to provide jagsUI-like summaries for NIMBLE output.
#
# Upstream licensing:
# - The BPAbook package is distributed under the GNU General Public License
#   version 3 (GPL-3). This derived/modified file is redistributed under GPL-3
#   as part of this repository. See the repository-level LICENSE file for the
#   full license text and terms.
#
# Modifications in this repository:
# - This file contains minor edits relative to the upstream BPAbook version.
#   Changes are documented in this repository's git history.
# ------------------------------------------------------------------------------
nimbleSummary <- function(samples, parameters = NULL) {

  # Handle WAIC
  WAIC <- NULL
  if (is.list(samples) && "samples" %in% names(samples)) {
    if ("WAIC" %in% names(samples)) WAIC <- samples$WAIC
    samples <- samples$samples
  }

  # Process into mcmc.list
  if (inherits(samples, "mcmc")) {
    samples <- coda::as.mcmc.list(samples)
  } else if (inherits(samples, "list")) {
    samples <- lapply(samples, coda::as.mcmc)
    samples <- coda::as.mcmc.list(samples)
  } else if (inherits(samples, "matrix")) {
    samples <- coda::as.mcmc(samples)
    samples <- coda::as.mcmc.list(samples)
  }
  if (!inherits(samples, "mcmc.list")) {
    stop("Unsupported input type", call. = FALSE)
  }

  # Reorder samples to match parameter vector
  if (!is.null(parameters)) {
    samples <- order_samples(samples, parameters)
  }

  out <- process_output(samples, DIC = FALSE, quiet = TRUE)
  out$WAIC <- WAIC
  out$samples <- samples
  nchains <- coda::nchain(samples)
  niter <- coda::niter(samples[[1]])
  out$mcmc.info <- list(n.chains = nchains, n.iter = niter,
                        n.samples = nchains * niter)
  class(out) <- c("nimbleSummary", "jagsUI")
  out
}

#------------------------------------------------------------------------------
# Get names of parameters from an mcmc.list
# If simplify = TRUE, also drop brackets/indices
param_names <- function(mcmc_list, simplify = FALSE) {
  out <- coda::varnames(mcmc_list)
  if (simplify) out <- strip_params(out, unique = TRUE)
  out
}

#------------------------------------------------------------------------------
# Match parameter name to scalar or array versions of parameter name
match_params <- function(params, params_raw) {
  unlist(lapply(params, function(x) {
    if (x %in% params_raw) return(x)
    if (!x %in% strip_params(params_raw)) return(NULL)
    params_raw[which_params(x, params_raw)]
  }))
}

#------------------------------------------------------------------------------
strip_params <- function(params_raw, unique = FALSE) {
  params_strip <- sapply(strsplit(params_raw, "[", fixed = TRUE), "[", 1)
  if (unique) return(unique(params_strip))
  params_strip
}

#------------------------------------------------------------------------------
# Identify which columns in mcmc.list object correspond to a given
# parameter name (useful for non-scalar parameters)
which_params <- function(param, params_raw) {
  params_strip <- strip_params(params_raw)
  if (!param %in% params_strip) {
    return(NULL)
  }
  which(params_strip == param)
}

#------------------------------------------------------------------------------
# Reorder output samples from coda to match input parameter order
order_samples <- function(samples, params) {
  tryCatch({
    matched <- match_params(params, param_names(samples))
    if ("deviance" %in% param_names(samples) & ! "deviance" %in% matched) {
      matched <- c(matched, "deviance")
    }
    samples[, matched, drop = FALSE]
  }, error = function(e) {
    message(paste0("Caught error re-ordering samples:\n", e, "\n"))
    samples
  })
}

#------------------------------------------------------------------------------
# Process output master function
# To generate backwards-compatible jagsUI output
process_output <- function(mcmc_list, coda_only = NULL, DIC, quiet = FALSE) {
  if (!quiet) {
    cat("Calculating statistics.......", "\n")
  }

  tryCatch({
    if (DIC == -999) stop("Throwing error for testing purposes", call. = FALSE)
    # Get the sims.list
    sims <- list(sims.list = sims_list(mcmc_list))
    # Calculate all stats
    stats <- calc_stats(mcmc_list, coda_only)
    # Convert them into stat arrays
    stats_list <- all_stat_arrays(stats, coda_only)
    # Get final summary table
    sum_list <- list(summary = stat_summary_table(stats, coda_only))
    # DIC stuff
    dic_list <- calc_DIC(mcmc_list, DIC)

    # Bind it all together
    if (!quiet) {
      cat("\nDone.", "\n")
    }
    c(sims, stats_list, dic_list, sum_list)
  }, error = function(e) {
    message(paste0("Processing output failed with this error:\n", e, "\n"))
    NULL
  })
}


#------------------------------------------------------------------------------
# Fill an array from vector using matching array indices
fill_array <- function(data_vector, indices) {
  out <- array(NA, dim = apply(indices, 2, max))
  out[indices] <- data_vector
  out
}


#------------------------------------------------------------------------------
# Extract the posterior of a parameter and organize it into an array
get_posterior_array <- function(parameter, samples) {

  tryCatch({
    #Subset output columns matching parameter
    col_inds <- which_params(parameter, param_names(samples))
    posterior_raw <- do.call(rbind, samples[, col_inds, drop = FALSE])

    #If parameter is scalar, return it now
    if (ncol(posterior_raw) == 1) {
      return(as.vector(posterior_raw))
    }

    #If parameter is array, get indices
    ind_raw <- get_inds(parameter, colnames(posterior_raw))
    ndraws <- nrow(posterior_raw)
    ind_array <- cbind(1:ndraws, ind_raw[rep(seq_len(nrow(ind_raw)), each = ndraws), ])

    #Create, fill, return output object
    fill_array(as.vector(posterior_raw), ind_array)
  }, error = function(e) {
    message(paste0("Caught error when creating sims.list array for '",
                   parameter, "':\n", e, "\n"))
    NA
  })
}


#------------------------------------------------------------------------------
# Get sims list
sims_list <- function(samples) {
  params <- param_names(samples)
  sapply(strip_params(params, unique = TRUE), get_posterior_array,
         samples, simplify = FALSE)
}


#------------------------------------------------------------------------------
# Extract stats for a parameter and organize into appropriately-sized array
get_stat_array <- function(parameter, stat, model_summary) {

  tryCatch({
    #Subset vector of stats for parameter
    row_ind <- which_params(parameter, rownames(model_summary))
    stat_vector <- model_summary[row_ind, stat]

    #If parameter is scalar, return it now
    if (length(stat_vector) == 1) return(stat_vector)

    #If parameter is array, get indices
    ind_array <- get_inds(parameter, names(stat_vector))

    #Create, fill, return output object
    fill_array(stat_vector, ind_array)
  }, error = function(e) {
    message(paste0("Caught error when creating stat array for '",
                   parameter, "':\n", e, "\n"))
    NA
  })
}


#------------------------------------------------------------------------------
# Compile all stats for all parameters into list of lists
all_stat_arrays <- function(summary_stats, coda_only) {

  stat_array_list <- function(stat, summary_stats) {
    params <- strip_params(rownames(summary_stats), unique = TRUE)
    sapply(params, function(x) {
      # If the parameter is in coda_only and the stat is not the mean, return NA
      if (x %in% coda_only && stat != "mean") return(NA)
      # Otherwise return the stat array for that parameter and stat
      get_stat_array(x, stat, summary_stats)
    }, simplify = FALSE)
  }
  # Do this for all stats
  out <- sapply(colnames(summary_stats), stat_array_list, summary_stats,
                simplify = FALSE)

  # Convert overlap0 to logical to match old jagsUI code
  out$overlap0 <- lapply(out$overlap0, function(x) x == 1)
  out
}


#------------------------------------------------------------------------------
# Convert stats into summary table in original jagsUI format
# For backwards compatibility
stat_summary_table <- function(stats, coda_only) {
  # Move overlap 0 and f to the end of the table
  stats <- stats[, c("mean", "sd", "q2.5", "q25", "q50", "q75", "q97.5",
                     "Rhat", "n.eff", "overlap0", "f"), drop = FALSE]
  # Rename the quantile columns
  colnames(stats)[3:7] <- c("2.5%", "25%", "50%", "75%", "97.5%")
  # Remove rows marked as coda_only
  keep_rows <- ! strip_params(rownames(stats)) %in% coda_only
  stats[keep_rows, , drop = FALSE]
}


#------------------------------------------------------------------------------
# Determine if 95% credible interval of parameter overlaps 0
overlap_0 <- function(lower, upper) {
  as.numeric(!(lower <= 0) == (upper < 0))
}

# Calculate proportion of posterior with same sign as mean
calc_f <- function(values, mn) {
  if (mn >= 0) return(mean(values >= 0, na.rm = TRUE))
  mean(values < 0, na.rm = TRUE)
}

calc_Rhat <- function(mcmc_list) {
  stopifnot(has_one_parameter(mcmc_list))
  if (length(mcmc_list) == 1) return(NA)
  out <- try(coda::gelman.diag(mcmc_list,
                               autoburnin = FALSE, multivariate = FALSE)$psrf[1])
  if (inherits(out, "try-error") || !is.finite(out)) out <- NA
  out
}

mcmc_to_mat <- function(mcmc_list) {
  stopifnot(has_one_parameter(mcmc_list))
  matrix(unlist(mcmc_list),
         nrow = coda::niter(mcmc_list), ncol = coda::nchain(mcmc_list))
}

# Based on R2WinBUGS code
calc_neff <- function(mcmc_list) {
  niter <- coda::niter(mcmc_list)
  nchain <- coda::nchain(mcmc_list)
  mcmc_mat <- mcmc_to_mat(mcmc_list)

  xdot <- apply(mcmc_mat, 2, mean, na.rm = TRUE)
  s2 <- apply(mcmc_mat, 2, stats::var, na.rm = TRUE)
  W <- mean(s2)

  # Non-degenerate case
  if (is.na(W)) {
    n_eff <- NA
  } else if ((W > 1.e-8) && (nchain > 1)) {
    B <- niter * stats::var(xdot)
    sig2hat <- ((niter - 1) * W + B) / niter
    n_eff <- round(nchain * niter * min(sig2hat / B, 1), 0)
  } else {
    # Degenerate case
    n_eff <- 1
  }
  n_eff
}

# Calculate series of statistics for one parameter
# Takes an mcmc.list as input
calc_param_stats <- function(mcmc_list, coda_only) {
  stopifnot(has_one_parameter(mcmc_list))
  values <- unlist(mcmc_list)
  stat_names <- c("mean", "sd", "q2.5", "q25", "q50", "q75", "q97.5",
                  "overlap0", "f", "Rhat", "n.eff")

  fallback <- sapply(stat_names, function(x) NA)
  if (any(is.infinite(values)) || all(is.na(values))) {
    return(fallback)
  }

  #Handle any unexpected errors during calculation
  tryCatch({
    # If the parameter is in codaOnly, return only the mean
    mn <- mean(values, na.rm = TRUE)
    if (coda_only) {
      fallback["mean"] <- mn
      return(fallback)
    }
    # Otherwise calculate all stats
    quants <- stats::quantile(values, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
    out <- c(mn,
             stats::sd(values, na.rm = TRUE),
             quants,
             overlap0 = overlap_0(quants[1], quants[5]),
             calc_f(values, mn),
             calc_Rhat(mcmc_list),
             calc_neff(mcmc_list))
    names(out) <- stat_names
    out
  }, error = function(e) {
    message(paste0("Caught error when calculating stats:\n", e, "\n"))
    fallback
  })
}


#------------------------------------------------------------------------------
# Calculate statistics for all parameters in posterior and organize into matrix
# Takes mcmc.list as input
calc_stats <- function(mcmc_list, coda_only = NULL) {
  params <- param_names(mcmc_list)
  coda_only <- strip_params(params) %in% coda_only

  out <- sapply(seq_along(params), function(i) {
    calc_param_stats(mcmc_list[, i], coda_only[i])
  })
  colnames(out) <- params
  t(out)
}


#------------------------------------------------------------------------------
# Calculate pD and DIC from deviance if it exists in output samples
calc_DIC <- function(samples, DIC) {
  if (!DIC || !("deviance" %in% param_names(samples))) {
    return(NULL)
  }

  dev <- mcmc_to_mat(samples[, "deviance"])
  if (any(is.na(dev)) || any(is.infinite(dev))) return(NULL)

  pd <- apply(dev, 2, FUN = function(x) stats::var(x) / 2)
  dic <- apply(dev, 2, mean) + pd

  c(pD = mean(pd), DIC = mean(dic))
}

#------------------------------------------------------------------------------
# Extract index values inside brackets from a non-scalar parameter
# param is the "base" name of the parameter and params_raw is a vector of
# strings that contain brackets
get_inds <- function(param, params_raw) {
  inds_raw <- sub(paste(param, "[", sep = ""), "", params_raw, fixed = TRUE)
  inds_raw <- sub("]", "", inds_raw, fixed = TRUE)
  inds_raw <- strsplit(inds_raw, ",", fixed = TRUE)
  inds <- as.integer(unlist(inds_raw))
  matrix(inds, byrow = TRUE, ncol = length(inds_raw[[1]]))
}


#------------------------------------------------------------------------------
# Check if mcmc.list has only one parameter (one column)
has_one_parameter <- function(mcmc_list) {
  coda::nvar(mcmc_list) == 1
}
