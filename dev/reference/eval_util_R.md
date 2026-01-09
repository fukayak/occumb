# Expected utility for regional species diversity assessments.

`eval_util_R()` evaluates the expected utility of a regional species
diversity assessment using Monte Carlo integration.

## Usage

``` r
eval_util_R(
  settings,
  fit = NULL,
  psi = NULL,
  theta = NULL,
  phi = NULL,
  N_rep = 1,
  cores = 1L
)
```

## Arguments

- settings:

  A data frame that specifies a set of conditions under which utility is
  evaluated. It must include columns named `J`, `K`, and `N`, which
  specify the number of sites, number of replicates per site, and
  sequencing depth per replicate, respectively. `J`, `K`, and `N` must
  be numeric vectors greater than 0. When `J` and `K` contain decimal
  values, they are discarded and treated as integers. Additional columns
  are ignored, but may be included.

- fit:

  An `occumbFit` object.

- psi:

  Sample values of the site occupancy probabilities of species stored in
  a matrix with sample \\\times\\ species dimensions or an array with
  sample \\\times\\ species \\\times\\ site dimensions.

- theta:

  Sample values of sequence capture probabilities of species stored in a
  matrix with sample \\\times\\ species dimensions or an array with
  sample \\\times\\ species \\\times\\ site dimensions.

- phi:

  Sample values of sequence relative dominance of species stored in a
  matrix with sample \\\times\\ species dimensions or an array with
  sample \\\times\\ species \\\times\\ site dimensions.

- N_rep:

  Controls the sample size for the Monte Carlo integration. The integral
  is evaluated using a total of `N_sample * N_rep` random samples, where
  `N_sample` is the maximum size of the MCMC sample in the `fit`
  argument and the parameter sample in the `psi`, `theta`, and `phi`
  arguments.

- cores:

  The number of cores to use for parallelization.

## Value

A data frame with a column named `Utility` in which the estimates of the
expected utility are stored. This is obtained by adding the `Utility`
column to the data frame provided in the `settings` argument.

## Details

The utility of a regional species diversity assessment can be defined as
the number of species expected to be detected in the region of interest
(Fukaya et al. 2022). `eval_util_R()` evaluates this utility for the
region modeled in the `occumbFit` object for the combination of `J`,
`K`, and `N` values specified in the `conditions` argument. Such
evaluations can be used to balance `J`, `K`, and `N` to maximize the
utility under a constant budget (possible combinations of `J`, `K`, and
`N` under a specified budget and cost values are easily obtained using
[`list_cond_R()`](https://fukayak.github.io/occumb/dev/reference/list_cond_R.md);
see the example below). It is also possible to examine how the utility
varies with different `J`, `K`, and `N` values without setting a budget
level, which may be useful in determining the satisfactory levels of
`J`, `K`, and `N` from a purely technical point of view.

The expected utility is defined as the expected value of the conditional
utility in the form: \$\$U(J, K, N \mid \boldsymbol{r}, \boldsymbol{u})
= \sum\_{i = 1}^{I}\left\\1 - \prod\_{j = 1}^{J}\prod\_{k =
1}^{K}\left(1 - \frac{u\_{ijk}r\_{ijk}}{\sum\_{m =
1}^{I}u\_{mjk}r\_{mjk}} \right)^N \right\\\$\$ where \\u\_{ijk}\\ is a
latent indicator variable representing the inclusion of the sequence of
species \\i\\ in replicate \\k\\ at site \\j\\, and \\r\_{ijk}\\ is a
latent variable that is proportional to the relative frequency of the
sequence of species \\i\\, conditional on its presence in replicate
\\k\\ at site \\j\\ (Fukaya et al. 2022). Expectations are taken with
respect to the posterior (or possibly prior) predictive distributions of
\\\boldsymbol{r} = \\r\_{ijk}\\\\ and \\\boldsymbol{u} = \\u\_{ijk}\\\\,
which are evaluated numerically using Monte Carlo integration. The
predictive distributions of \\\boldsymbol{r}\\ and \\\boldsymbol{u}\\
depend on the model parameters \\\psi\\, \\\theta\\, and \\\phi\\
values. Their posterior (or prior) distribution is specified by
supplying an `occumbFit` object containing their posterior samples via
the `fit` argument, or by supplying a matrix or array of posterior (or
prior) samples of parameter values via the `psi`, `theta`, and `phi`
arguments. Higher approximation accuracy can be obtained by increasing
the value of `N_rep`.

The `eval_util_R()` function can be executed by supplying the `fit`
argument without specifying the `psi`, `theta`, and `phi` arguments, by
supplying the three `psi`, `theta`, and `phi` arguments without the
`fit` argument, or by supplying the `fit` argument and any or all of the
`psi`, `theta`, and `phi` arguments. If the `psi`, `theta`, or `phi`
arguments are specified in addition to the `fit`, the parameter values
given in these arguments are preferentially used to evaluate the
expected utility. If the sample sizes differed among parameters,
parameters with smaller sample sizes are resampled with replacements to
align the sample sizes across parameters.

The expected utility is evaluated assuming homogeneity of replicates, in
the sense that \\\theta\\ and \\\phi\\, the model parameters associated
with the species detection process, are constant across replicates
within a site. For this reason, `eval_util_R()` does not accept
replicate-specific \\\theta\\ and \\\phi\\. If the `occumbFit` object
supplied in the `fit` argument has a replicate-specific parameter, the
parameter samples to be used in the utility evaluation must be provided
explicitly via the `theta` or `phi` arguments.

If the parameters are modeled as a function of site covariates in the
`fit` object, or if the `psi`, `theta`, and/or `phi` arguments have site
dimensions, the expected utility is evaluated to account for the site
heterogeneity of the parameters. To incorporate site heterogeneity, the
parameter values for each `J` site are determined by selecting
site-specific parameter values in the `fit`, or those supplied in `psi`,
`theta`, and `phi` via random sampling with replacement. Thus, expected
utility is evaluated by assuming a set of supplied parameter values as a
statistical population of site-specific parameters.

The Monte Carlo integration is executed in parallel on multiple CPU
cores, where the `cores` argument controls the degree of
parallelization.

## References

K. Fukaya, N. I. Kondo, S. S. Matsuzaki and T. Kadoya (2022)
Multispecies site occupancy modelling and study design for spatially
replicated environmental DNA metabarcoding. *Methods in Ecology and
Evolution* **13**:183–193.
[doi:10.1111/2041-210X.13732](https://doi.org/10.1111/2041-210X.13732)

## Examples

``` r
# \donttest{
set.seed(1)

# Generate a random dataset (20 species * 2 sites * 2 reps)
I <- 20 # Number of species
J <- 2  # Number of sites
K <- 2  # Number of replicates
data <- occumbData(
    y = array(sample.int(I * J * K), dim = c(I, J, K)))

# Fitting a null model
fit <- occumb(data = data)
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
#>    Unobserved stochastic nodes: 229
#>    Total graph size: 673
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

## Estimate expected utility
# Arbitrary J, K, and N values
(util1 <- eval_util_R(expand.grid(J = 1:3, K = 1:3, N = c(1E3, 1E4, 1E5)),
                      fit))
#>    J K     N  Utility
#> 1  1 1 1e+03 19.80640
#> 2  2 1 1e+03 19.97425
#> 3  3 1 1e+03 19.99291
#> 4  1 2 1e+03 19.88613
#> 5  2 2 1e+03 19.97950
#> 6  3 2 1e+03 19.99450
#> 7  1 3 1e+03 19.88699
#> 8  2 3 1e+03 19.98400
#> 9  3 3 1e+03 19.99450
#> 10 1 1 1e+04 19.84006
#> 11 2 1 1e+04 19.97569
#> 12 3 1 1e+04 19.99100
#> 13 1 2 1e+04 19.89600
#> 14 2 2 1e+04 19.97875
#> 15 3 2 1e+04 19.99200
#> 16 1 3 1e+04 19.89425
#> 17 2 3 1e+04 19.98425
#> 18 3 3 1e+04 19.99150
#> 19 1 1 1e+05 19.84624
#> 20 2 1 1e+05 19.97525
#> 21 3 1 1e+05 19.99075
#> 22 1 2 1e+05 19.88575
#> 23 2 2 1e+05 19.97850
#> 24 3 2 1e+05 19.99250
#> 25 1 3 1e+05 19.88675
#> 26 2 3 1e+05 19.97975
#> 27 3 3 1e+05 19.99325

# J, K, and N values under specified budget and cost
(util2 <- eval_util_R(list_cond_R(budget = 50000,
                                  lambda1 = 0.01,
                                  lambda2 = 5000,
                                  lambda3 = 5000),
                      fit))
#>    budget lambda1 lambda2 lambda3 J K          N  Utility
#> 1   50000    0.01    5000    5000 1 1 4000000.00 19.83575
#> 2   50000    0.01    5000    5000 2 1 1500000.00 19.97650
#> 3   50000    0.01    5000    5000 3 1  666666.67 19.99425
#> 4   50000    0.01    5000    5000 4 1  250000.00 19.99675
#> 5   50000    0.01    5000    5000 1 2 1750000.00 19.89075
#> 6   50000    0.01    5000    5000 2 2  500000.00 19.98075
#> 7   50000    0.01    5000    5000 3 2   83333.33 19.99250
#> 8   50000    0.01    5000    5000 1 3 1000000.00 19.89250
#> 9   50000    0.01    5000    5000 2 3  166666.67 19.98100
#> 10  50000    0.01    5000    5000 1 4  625000.00 19.89600
#> 11  50000    0.01    5000    5000 1 5  400000.00 19.88925
#> 12  50000    0.01    5000    5000 1 6  250000.00 19.90025
#> 13  50000    0.01    5000    5000 1 7  142857.14 19.89300
#> 14  50000    0.01    5000    5000 1 8   62500.00 19.88900

# K values restricted
(util3 <- eval_util_R(list_cond_R(budget = 50000,
                                  lambda1 = 0.01,
                                  lambda2 = 5000,
                                  lambda3 = 5000,
                                  K = 1:5),
                      fit))
#>    budget lambda1 lambda2 lambda3 J K          N  Utility
#> 1   50000    0.01    5000    5000 1 1 4000000.00 19.85850
#> 2   50000    0.01    5000    5000 2 1 1500000.00 19.97550
#> 3   50000    0.01    5000    5000 3 1  666666.67 19.99100
#> 4   50000    0.01    5000    5000 4 1  250000.00 19.99500
#> 5   50000    0.01    5000    5000 1 2 1750000.00 19.90200
#> 6   50000    0.01    5000    5000 2 2  500000.00 19.97850
#> 7   50000    0.01    5000    5000 3 2   83333.33 19.99400
#> 8   50000    0.01    5000    5000 1 3 1000000.00 19.88375
#> 9   50000    0.01    5000    5000 2 3  166666.67 19.98100
#> 10  50000    0.01    5000    5000 1 4  625000.00 19.89225
#> 11  50000    0.01    5000    5000 1 5  400000.00 19.89425

# J and K values restricted
(util4 <- eval_util_R(list_cond_R(budget = 50000,
                                  lambda1 = 0.01,
                                  lambda2 = 5000,
                                  lambda3 = 5000,
                                  J = 1:3, K = 1:5),
                      fit))
#>    budget lambda1 lambda2 lambda3 J K          N  Utility
#> 1   50000    0.01    5000    5000 1 1 4000000.00 19.84800
#> 2   50000    0.01    5000    5000 2 1 1500000.00 19.97450
#> 3   50000    0.01    5000    5000 3 1  666666.67 19.99100
#> 4   50000    0.01    5000    5000 1 2 1750000.00 19.88450
#> 5   50000    0.01    5000    5000 2 2  500000.00 19.97850
#> 6   50000    0.01    5000    5000 3 2   83333.33 19.99150
#> 7   50000    0.01    5000    5000 1 3 1000000.00 19.89500
#> 8   50000    0.01    5000    5000 2 3  166666.67 19.97900
#> 9   50000    0.01    5000    5000 1 4  625000.00 19.89375
#> 10  50000    0.01    5000    5000 1 5  400000.00 19.89575

# theta and phi values supplied
(util5 <- eval_util_R(list_cond_R(budget = 50000,
                                  lambda1 = 0.01,
                                  lambda2 = 5000,
                                  lambda3 = 5000,
                                  J = 1:3, K = 1:5),
                      fit,
                      theta = array(0.5, dim = c(4000, I, J)),
                      phi = array(1, dim = c(4000, I, J))))
#>    budget lambda1 lambda2 lambda3 J K          N  Utility
#> 1   50000    0.01    5000    5000 1 1 4000000.00  9.97825
#> 2   50000    0.01    5000    5000 2 1 1500000.00 14.90850
#> 3   50000    0.01    5000    5000 3 1  666666.67 17.42767
#> 4   50000    0.01    5000    5000 1 2 1750000.00 14.88425
#> 5   50000    0.01    5000    5000 2 2  500000.00 18.71124
#> 6   50000    0.01    5000    5000 3 2   83333.33 19.66727
#> 7   50000    0.01    5000    5000 1 3 1000000.00 17.42300
#> 8   50000    0.01    5000    5000 2 3  166666.67 19.65025
#> 9   50000    0.01    5000    5000 1 4  625000.00 18.66550
#> 10  50000    0.01    5000    5000 1 5  400000.00 19.25875

# psi, theta, and phi values, but no fit object supplied
(util6 <- eval_util_R(list_cond_R(budget = 50000,
                                  lambda1 = 0.01,
                                  lambda2 = 5000,
                                  lambda3 = 5000,
                                  J = 1:3, K = 1:5),
                      fit = NULL,
                      psi = array(0.9, dim = c(4000, I, J)),
                      theta = array(0.9, dim = c(4000, I, J)),
                      phi = array(1, dim = c(4000, I, J))))
#>    budget lambda1 lambda2 lambda3 J K          N  Utility
#> 1   50000    0.01    5000    5000 1 1 4000000.00 16.18223
#> 2   50000    0.01    5000    5000 2 1 1500000.00 19.28467
#> 3   50000    0.01    5000    5000 3 1  666666.67 19.86450
#> 4   50000    0.01    5000    5000 1 2 1750000.00 17.79325
#> 5   50000    0.01    5000    5000 2 2  500000.00 19.76721
#> 6   50000    0.01    5000    5000 3 2   83333.33 19.97950
#> 7   50000    0.01    5000    5000 1 3 1000000.00 17.99450
#> 8   50000    0.01    5000    5000 2 3  166666.67 19.80000
#> 9   50000    0.01    5000    5000 1 4  625000.00 17.97675
#> 10  50000    0.01    5000    5000 1 5  400000.00 17.98575
# }
```
