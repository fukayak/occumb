# occumb: Site Occupancy Modeling For Environmental DNA Metabarcoding

<!-- badges: start -->
[![R-CMD-check](https://github.com/fukayak/occumb/workflows/R-CMD-check/badge.svg)](https://github.com/fukayak/occumb/actions)
<!-- badges: end -->

occumb is an R package that aims to apply [the multispecies site occupancy modeling for environmental DNA (eDNA) metabarcoding](https://doi.org/10.1111/2041-210X.13732) easily.

This package allows the users to fit the model with a fully Bayesian approach, using the conventional formulas in R. It enables the analysis of the detectability of species in different stages of the workflow of eDNA metabarcoding and inference of species site occupancy while accounting for false negatives. It also provides the functionality for a model-based inference to assist the optimization of the study design.

See the [original paper](https://doi.org/10.1111/2041-210X.13732) for details of the model and inference.

## Important note
The current programs are an alpha version released to present the concept of package development. In other words, please be aware of the following points when you want to try this package:

- Some of the planned features have not yet been implemented.
- The functions may not have been thoroughly tested.
- Documentation may be insufficient.
- Specifications may change significantly in the future.

## Installation

You can install the latest version of occumb from the GitHub repository using the `install_github()` function in remotes (or devtools) package:

``` r
remotes::install_github("fukayak/occumb", ref = "main")
```

## Example

The current version (v0.4.x) provides the following main functions:

- `occumbData()` builds a dataset object used for the multispecies site occupancy modeling.
- `occumb()` fits a model, possibly with species, site, and replicate covariates.
- `gof()` computes statistics for goodness-of-fit assessment of the model.
- `loglik()` extracts the pointwise log-likelihood of the model.
- `eval_util_L()` and `eval_util_R()` predicts expected utility (the expected number of detected species) at the local and regional scale, respectively.

``` r
library(occumb)

# Generate the smallest random dataset (2 species * 2 sites * 2 reps)
I <- 2 # Number of species
J <- 2 # Number of sites
K <- 2 # Number of replicates
data <- occumbData(
    y = array(sample.int(I * J * K), dim = c(I, J, K)),
    spec_cov = list(cov1 = rnorm(I)),
    site_cov = list(cov2 = rnorm(J),
                    cov3 = factor(1:J)),
    repl_cov = list(cov4 = matrix(rnorm(J * K), J, K)))

# Fitting a null model (includes only species-specific intercepts)
res0 <- occumb(data = data)

# occumb() fits the model using jags() in the jagsUI package
# You can thus get the same output as seen in the jagsUI results
res0
summary(res0)
plot(res0)

# Add species-specific effects of site covariates in occupancy probabilities
res1 <- occumb(formula_psi = ~ cov2, data = data)        # Continuous covariate
res2 <- occumb(formula_psi = ~ cov3, data = data)        # Categorical covariate
res3 <- occumb(formula_psi = ~ cov2 * cov3, data = data) # Interaction

# Add species covariate in the three parameters
# Note that species covariates are modeled as common effects
res4 <- occumb(formula_phi_shared = ~ cov1, data = data)   # phi
res5 <- occumb(formula_theta_shared = ~ cov1, data = data) # theta
res6 <- occumb(formula_psi_shared = ~ cov1, data = data)   # psi

# Add replicate covariates
# Note that replicate covariates can only be specified for theta and phi
res7 <- occumb(formula_phi = ~ cov4, data = data)   # phi
res8 <- occumb(formula_theta = ~ cov4, data = data) # theta

# Specify the prior distribution and MCMC settings explicitly
res9 <- occumb(data = data, prior_prec = 1E-2, prior_ulim = 1E2,
               n.chains = 1, n.burnin = 1000, n.thin = 1, n.iter = 2000)
res10 <- occumb(data = data, parallel = TRUE) # Run MCMC in parallel

# Calculate goodness-of-fit statistics
gof_res0 <- gof(res0)

# Extract pointwise log-likelihood
ll_res0 <- loglik(res0)

# Estimate expected utility for local species diversity assessment
util_L <- eval_util_L(list_cond_L(budget = 1E5, lambda1 = 0.01, lambda2 = 5000, res0),
                      res0)

# Estimate expected utility for regional species diversity assessment
util_R <- eval_util_R(list_cond_R(budget = 50000,
                                  lambda1 = 0.01,
                                  lambda2 = 5000,
                                  lambda3 = 5000),
                      res0)
```

See the documentation of each function for more details.

