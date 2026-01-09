# Model-fitting function.

`occumb()` fits the community site-occupancy model for eDNA
metabarcoding (Fukaya et al. 2022) and returns a model-fit object
containing posterior samples.

## Usage

``` r
occumb(
  formula_phi = ~1,
  formula_theta = ~1,
  formula_psi = ~1,
  formula_phi_shared = ~1,
  formula_theta_shared = ~1,
  formula_psi_shared = ~1,
  prior_prec = 1e-04,
  prior_ulim = 10000,
  data,
  n.chains = 4,
  n.adapt = NULL,
  n.burnin = 10000,
  n.thin = 10,
  n.iter = 20000,
  parallel = FALSE,
  ...
)
```

## Arguments

- formula_phi:

  A right-hand side formula describing species-specific effects of
  sequence relative dominance (\\\phi\\).

- formula_theta:

  A right-hand side formula describing species-specific effects of
  sequence capture probability (\\\theta\\).

- formula_psi:

  A right-hand side formula describing species-specific effects of
  occupancy probability (\\\psi\\).

- formula_phi_shared:

  A right-hand side formula describing effects of sequence relative
  dominance (\\\phi\\) that are common across species. The intercept
  term is ignored (see Details).

- formula_theta_shared:

  A right-hand side formula describing effects of sequence capture
  probability (\\\theta\\) that are common across species. The intercept
  term is ignored (see Details).

- formula_psi_shared:

  A right-hand side formula describing effects of occupancy probability
  (\\\psi\\) that are common across species. The intercept term is
  ignored (see Details).

- prior_prec:

  Precision of normal prior distribution for the community-level average
  of species-specific parameters and effects common across species.

- prior_ulim:

  Upper limit of uniform prior distribution for the standard deviation
  of species-specific parameters.

- data:

  A dataset supplied as an
  [`occumbData`](https://fukayak.github.io/occumb/reference/occumbData.md)
  class object.

- n.chains:

  Number of Markov chains to run.

- n.adapt:

  Number of iterations to run in the JAGS adaptive phase.

- n.burnin:

  Number of iterations at the beginning of the chain to discard.

- n.thin:

  Thinning rate. Must be a positive integer.

- n.iter:

  Total number of iterations per chain (including burn-in).

- parallel:

  If TRUE, run MCMC chains in parallel on multiple CPU cores.

- ...:

  Additional arguments passed to
  [`jags()`](https://kenkellner.com/jagsUI/reference/jags.html)
  function.

## Value

An S4 object of the `occumbFit` class containing the results of the
model fitting and the supplied dataset.

## Details

`occumb()` allows the fitting of a range of community site occupancy
models, including covariates at different levels of the data generation
process. The most general form of the model can be written as follows
(the notation follows that of the original article; see Fukaya et al.
(2022) or [the package
vignette](https://fukayak.github.io/occumb/articles/model_specification.html)).

Sequence read counts: \$\$(y\_{1jk}, ..., y\_{Ijk}) \sim
\textrm{Multinomial}((\pi\_{1jk}, ..., \pi\_{Ijk}), N\_{jk}),\$\$
\$\$\pi\_{ijk} = \frac{u\_{ijk}r\_{ijk}}{\sum_m u\_{mjk}r\_{mjk}},\$\$

Relative frequency of species sequences: \$\$r\_{ijk} \sim
\textrm{Gamma}(\phi\_{ijk}, 1),\$\$

Capture of species sequences: \$\$u\_{ijk} \sim
\textrm{Bernoulli}(z\_{ij}\theta\_{ijk}),\$\$

Site occupancy of species: \$\$z\_{ij} \sim
\textrm{Bernoulli}(\psi\_{ij}),\$\$ where the variations of \\\phi\\,
\\\theta\\, and \\\psi\\ are modeled by specifying model formulas in
`formula_phi`, `formula_theta`, `formula_psi`, `formula_phi_shared`,
`formula_theta_shared`, and `formula_psi_shared.` Each parameter may
have species-specific effects and effects that are common across
species, where the former is specified by `formula_phi`,
`formula_theta`, and `formula_psi`, whereas `formula_phi_shared`,
`formula_theta_shared`, and `formula_psi_shared` specify the latter. As
species-specific intercepts are specified by default, the intercept
terms in `formula_phi_shared`, `formula_theta_shared`, and
`formula_psi_shared` are always ignored. Covariate terms must be found
in the names of the list elements stored in the `spec_cov`, `site_cov`,
or `repl_cov` slots in the dataset object provided with the `data`
argument. Covariates are modeled using the log link function for
\\\phi\\ and logit link function for \\\theta\\ and \\\psi.\\

The two arguments, `prior_prec` and `prior_ulim`, control the prior
distribution of parameters. For the community-level average of
species-specific effects and effects common across species, a normal
prior distribution with a mean of 0 and precision (i.e., the inverse of
the variance) `prior_prec` is specified. For the standard deviation of
species-specific effects, a uniform prior distribution with a lower
limit of zero and an upper limit of `prior_ulim` is specified. For the
correlation coefficient of species-specific effects, a uniform prior
distribution in the range of \\-\\1 to 1 is specified by default.

See [the package
vignette](https://fukayak.github.io/occumb/articles/model_specification.html)
for details on the model specifications in `occumb()`.

The `data` argument requires a dataset object to be generated using
`ocumbData()`; see the document of
[`occumbData()`](https://fukayak.github.io/occumb/reference/occumbData.md).

The model is fit using the
[`jags()`](https://kenkellner.com/jagsUI/reference/jags.html) function
of the [jagsUI](https://cran.r-project.org/package=jagsUI) package,
where Markov chain Monte Carlo (MCMC) methods are used to obtain
posterior samples of the parameters and latent variables. Arguments
`n.chains`, `n.adapt`, `n.burnin`, `n.thin`, `n.iter`, and `parallel`
are passed on to arguments of the same name in the
[`jags()`](https://kenkellner.com/jagsUI/reference/jags.html) function.
See the document of
[jagsUI](https://cran.r-project.org/package=jagsUI)'s
[`jags()`](https://kenkellner.com/jagsUI/reference/jags.html) function
for details. A set of random initial values is used to perform an MCMC
run.

## References

K. Fukaya, N. I. Kondo, S. S. Matsuzaki and T. Kadoya (2022)
Multispecies site occupancy modelling and study design for spatially
replicated environmental DNA metabarcoding. *Methods in Ecology and
Evolution* **13**:183–193.
[doi:10.1111/2041-210X.13732](https://doi.org/10.1111/2041-210X.13732)

## Examples

``` r
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

# \donttest{
# Fitting a null model (includes only species-specific intercepts)
res0 <- occumb(data = data)
#> 
#> Processing function input....... 
#> 
#> Done. 
#>  
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 31
#>    Total graph size: 133
#> 
#> Initializing model
#> 
#> Adaptive phase..... 
#> Adaptive phase complete 
#>  
#> 
#>  Burn-in phase, 10000 iterations x 4 chains 
#>  
#> 
#> Sampling from joint posterior, 10000 iterations x 4 chains 
#>  
#> 
#> Calculating statistics....... 
#> 
#> Done. 

# Add species-specific effects of site covariates in occupancy probabilities
res1 <- occumb(formula_psi = ~ cov2, data = data)        # Continuous covariate
#> 
#> Processing function input....... 
#> 
#> Done. 
#>  
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 36
#>    Total graph size: 156
#> 
#> Initializing model
#> 
#> Adaptive phase..... 
#> Adaptive phase complete 
#>  
#> 
#>  Burn-in phase, 10000 iterations x 4 chains 
#>  
#> 
#> Sampling from joint posterior, 10000 iterations x 4 chains 
#>  
#> 
#> Calculating statistics....... 
#> 
#> Done. 
res2 <- occumb(formula_psi = ~ cov3, data = data)        # Categorical covariate
#> 
#> Processing function input....... 
#> 
#> Done. 
#>  
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 36
#>    Total graph size: 156
#> 
#> Initializing model
#> 
#> Adaptive phase..... 
#> Adaptive phase complete 
#>  
#> 
#>  Burn-in phase, 10000 iterations x 4 chains 
#>  
#> 
#> Sampling from joint posterior, 10000 iterations x 4 chains 
#>  
#> 
#> Calculating statistics....... 
#> 
#> Done. 
res3 <- occumb(formula_psi = ~ cov2 * cov3, data = data) # Interaction
#> 
#> Processing function input....... 
#> 
#> Done. 
#>  
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 49
#>    Total graph size: 190
#> 
#> Initializing model
#> 
#> Adaptive phase..... 
#> Adaptive phase complete 
#>  
#> 
#>  Burn-in phase, 10000 iterations x 4 chains 
#>  
#> 
#> Sampling from joint posterior, 10000 iterations x 4 chains 
#>  
#> 
#> Calculating statistics....... 
#> 
#> Done. 
# Add species covariate in the three parameters
# Note that species covariates are modeled as common effects
res4 <- occumb(formula_phi_shared = ~ cov1, data = data)   # phi
#> 
#> Processing function input....... 
#> 
#> Done. 
#>  
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 32
#>    Total graph size: 141
#> 
#> Initializing model
#> 
#> Adaptive phase..... 
#> Adaptive phase complete 
#>  
#> 
#>  Burn-in phase, 10000 iterations x 4 chains 
#>  
#> 
#> Sampling from joint posterior, 10000 iterations x 4 chains 
#>  
#> 
#> Calculating statistics....... 
#> 
#> Done. 
res5 <- occumb(formula_theta_shared = ~ cov1, data = data) # theta
#> 
#> Processing function input....... 
#> 
#> Done. 
#>  
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 32
#>    Total graph size: 141
#> 
#> Initializing model
#> 
#> Adaptive phase..... 
#> Adaptive phase complete 
#>  
#> 
#>  Burn-in phase, 10000 iterations x 4 chains 
#>  
#> 
#> Sampling from joint posterior, 10000 iterations x 4 chains 
#>  
#> 
#> Calculating statistics....... 
#> 
#> Done. 
res6 <- occumb(formula_psi_shared = ~ cov1, data = data)   # psi
#> 
#> Processing function input....... 
#> 
#> Done. 
#>  
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 32
#>    Total graph size: 141
#> 
#> Initializing model
#> 
#> Adaptive phase..... 
#> Adaptive phase complete 
#>  
#> 
#>  Burn-in phase, 10000 iterations x 4 chains 
#>  
#> 
#> Sampling from joint posterior, 10000 iterations x 4 chains 
#>  
#> 
#> Calculating statistics....... 
#> 
#> Done. 
# Add replicate covariates
# Note that replicate covariates can only be specified for theta and phi
res7 <- occumb(formula_phi = ~ cov4, data = data)   # phi
#> 
#> Processing function input....... 
#> 
#> Done. 
#>  
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 36
#>    Total graph size: 170
#> 
#> Initializing model
#> 
#> Adaptive phase..... 
#> Adaptive phase complete 
#>  
#> 
#>  Burn-in phase, 10000 iterations x 4 chains 
#>  
#> 
#> Sampling from joint posterior, 10000 iterations x 4 chains 
#>  
#> 
#> Calculating statistics....... 
#> 
#> Done. 
res8 <- occumb(formula_theta = ~ cov4, data = data) # theta
#> 
#> Processing function input....... 
#> 
#> Done. 
#>  
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 36
#>    Total graph size: 174
#> 
#> Initializing model
#> 
#> Adaptive phase..... 
#> Adaptive phase complete 
#>  
#> 
#>  Burn-in phase, 10000 iterations x 4 chains 
#>  
#> 
#> Sampling from joint posterior, 10000 iterations x 4 chains 
#>  
#> 
#> Calculating statistics....... 
#> 
#> Done. 
# Specify the prior distribution and MCMC settings explicitly
res9 <- occumb(data = data, prior_prec = 1E-2, prior_ulim = 1E2,
               n.chains = 1, n.burnin = 1000, n.thin = 1, n.iter = 2000)
#> 
#> Processing function input....... 
#> 
#> Done. 
#>  
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 31
#>    Total graph size: 133
#> 
#> Initializing model
#> 
#> Adaptive phase..... 
#> Adaptive phase complete 
#>  
#> 
#>  Burn-in phase, 1000 iterations x 1 chains 
#>  
#> 
#> Sampling from joint posterior, 1000 iterations x 1 chains 
#>  
#> 
#> Calculating statistics....... 
#> 
#> Done. 
res10 <- occumb(data = data, parallel = TRUE, n.cores = 2) # Run MCMC in parallel
#> 
#> Processing function input....... 
#> 
#> Done. 
#>  
#> Beginning parallel processing using 2 cores. Console output will be suppressed.
#> 
#> Parallel processing completed.
#> 
#> Calculating statistics....... 
#> 
#> Done. 
# }
```
