# Expected utility for local species diversity assessments.

`eval_util_L()` evaluates the expected utility of a local species
diversity assessment by using Monte Carlo integration.

## Usage

``` r
eval_util_L(
  settings,
  fit = NULL,
  z = NULL,
  theta = NULL,
  phi = NULL,
  N_rep = 1,
  cores = 1L
)
```

## Arguments

- settings:

  A data frame that specifies a set of conditions under which utility is
  evaluated. It must include columns named `K` and `N`, which specify
  the number of replicates per site and the sequencing depth per
  replicate, respectively. `K` and `N` must be numeric vectors greater
  than 0. When `K` contains a decimal value, it is discarded and treated
  as an integer. Additional columns are ignored, but may be included.

- fit:

  An `occumbFit` object.

- z:

  Sample values of site occupancy status of species stored in an array
  with sample \\\times\\ species \\\times\\ site dimensions.

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
  is evaluated using `N_sample * N_rep` random samples, where `N_sample`
  is the maximum size of the MCMC sample in the `fit` argument and the
  parameter sample in the `z`, `theta`, and `phi` arguments.

- cores:

  The number of cores to use for parallelization.

## Value

A data frame with a column named `Utility` in which the estimates of the
expected utility are stored. This is obtained by adding the `Utility`
column to the data frame provided in the `settings` argument.

## Details

The utility of local species diversity assessment for a given set of
sites can be defined as the expected number of detected species per site
(Fukaya et al. 2022). `eval_util_L()` evaluates this utility for
arbitrary sets of sites that can potentially have different values for
site occupancy status of species, \\z\\, sequence capture probabilities
of species, \\\theta\\, and sequence relative dominance of species,
\\\phi\\, for the combination of `K` and `N` values specified in the
`conditions` argument. Such evaluations can be used to balance `K` and
`N` to maximize the utility under a constant budget (possible
combinations of `K` and `N` under a specified budget and cost values are
easily obtained using
[`list_cond_L()`](https://fukayak.github.io/occumb/reference/list_cond_L.md);
see the example below). It is also possible to examine how the utility
varies with different `K` and `N` values without setting a budget level,
which may be useful for determining a satisfactory level of `K` and `N`
from a purely technical point of view. The expected utility is defined
as the expected value of the conditional utility in the form: \$\$U(K, N
\mid \boldsymbol{r}, \boldsymbol{u}) = \frac{1}{J}\sum\_{j =
1}^{J}\sum\_{i = 1}^{I}\left\\1 - \prod\_{k = 1}^{K}\left(1 -
\frac{u\_{ijk}r\_{ijk}}{\sum\_{m = 1}^{I}u\_{mjk}r\_{mjk}} \right)^N
\right\\\$\$ where \\u\_{ijk}\\ is a latent indicator variable
representing the inclusion of the sequence of species \\i\\ in replicate
\\k\\ at site \\j\\, and \\r\_{ijk}\\ is a latent variable that is
proportional to the relative frequency of the sequence of species \\i\\,
conditional on its presence in replicate \\k\\ at site \\j\\ (Fukaya et
al. 2022). Expectations are taken with respect to the posterior (or
possibly prior) predictive distributions of \\\boldsymbol{r} =
\\r\_{ijk}\\\\ and \\\boldsymbol{u} = \\u\_{ijk}\\\\, which are
evaluated numerically using Monte Carlo integration. The predictive
distributions of \\\boldsymbol{r}\\ and \\\boldsymbol{u}\\ depend on the
model parameters \\z\\, \\\theta\\, and \\\phi\\ values. Their posterior
(or prior) distribution is specified by supplying an `occumbFit` object
containing their posterior samples via the `fit` argument, or by
supplying a matrix or array of posterior (or prior) samples of parameter
values via the `z`, `theta`, and `phi` arguments. Higher approximation
accuracy can be obtained by increasing the value of `N_rep`.

The `eval_util_L()` function can be executed by supplying the `fit`
argument without specifying the `z`, `theta`, and `phi` arguments, by
supplying the three `z`, `theta`, and `phi` arguments without the `fit`
argument, or by supplying the `fit` argument and any or all of the `z`,
`theta`, and `phi` arguments. If `z`, `theta`, or `phi` arguments are
specified in addition to the `fit`, the parameter values given in these
arguments are used preferentially to evaluate the expected utility. If
the sample sizes differ among parameters, parameters with smaller sample
sizes are resampled with replacements to align the sample sizes across
parameters.

The expected utility is evaluated assuming homogeneity of replicates, in
the sense that \\\theta\\ and \\\phi\\, the model parameters associated
with the species detection process, are constant across replicates
within a site. For this reason, `eval_util_L()` does not accept
replicate-specific \\\theta\\ and \\\phi\\. If the `occumbFit` object
supplied in the `fit` argument has a replicate-specific parameter, the
parameter samples to be used in the utility evaluation must be provided
explicitly via the `theta` or `phi` arguments.

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
# Arbitrary K and N values
(util1 <- eval_util_L(expand.grid(K = 1:3, N = c(1E3, 1E4, 1E5)),
                      fit))
#>   K     N  Utility
#> 1 1 1e+03 19.92034
#> 2 2 1e+03 19.99545
#> 3 3 1e+03 19.99875
#> 4 1 1e+04 19.95207
#> 5 2 1e+04 19.99525
#> 6 3 1e+04 19.99888
#> 7 1 1e+05 19.95416
#> 8 2 1e+05 19.99475
#> 9 3 1e+05 19.99875

# K and N values under specified budget and cost
(util2 <- eval_util_L(list_cond_L(budget = 1E5,
                                  lambda1 = 0.01,
                                  lambda2 = 5000,
                                  fit),
                      fit))
#>   budget lambda1 lambda2 K          N  Utility
#> 1  1e+05    0.01    5000 1 4500000.00 19.95200
#> 2  1e+05    0.01    5000 2 2000000.00 19.99588
#> 3  1e+05    0.01    5000 3 1166666.67 19.99912
#> 4  1e+05    0.01    5000 4  750000.00 19.99975
#> 5  1e+05    0.01    5000 5  500000.00 19.99963
#> 6  1e+05    0.01    5000 6  333333.33 20.00000
#> 7  1e+05    0.01    5000 7  214285.71 20.00000
#> 8  1e+05    0.01    5000 8  125000.00 20.00000
#> 9  1e+05    0.01    5000 9   55555.56 20.00000

# K values restricted
(util3 <- eval_util_L(list_cond_L(budget = 1E5,
                                  lambda1 = 0.01,
                                  lambda2 = 5000,
                                  fit,
                                  K = 1:5),
                      fit))
#>   budget lambda1 lambda2 K       N  Utility
#> 1  1e+05    0.01    5000 1 4500000 19.95175
#> 2  1e+05    0.01    5000 2 2000000 19.99600
#> 3  1e+05    0.01    5000 3 1166667 19.99850
#> 4  1e+05    0.01    5000 4  750000 19.99963
#> 5  1e+05    0.01    5000 5  500000 19.99963

# theta and phi values supplied
(util4 <- eval_util_L(list_cond_L(budget = 1E5,
                                  lambda1 = 0.01,
                                  lambda2 = 5000,
                                  fit,
                                  K = 1:5),
                      fit,
                      theta = array(0.5, dim = c(4000, I, J)),
                      phi = array(1, dim = c(4000, I, J))))
#>   budget lambda1 lambda2 K       N  Utility
#> 1  1e+05    0.01    5000 1 4500000  9.96600
#> 2  1e+05    0.01    5000 2 2000000 15.00973
#> 3  1e+05    0.01    5000 3 1166667 17.49592
#> 4  1e+05    0.01    5000 4  750000 18.73355
#> 5  1e+05    0.01    5000 5  500000 19.36387

# z, theta, and phi values, but no fit object supplied
(util5 <- eval_util_L(list_cond_L(budget = 1E5,
                                  lambda1 = 0.01,
                                  lambda2 = 5000,
                                  fit,
                                  K = 1:5),
                      fit = NULL,
                      z = array(1, dim = c(4000, I, J)),
                      theta = array(0.5, dim = c(4000, I, J)),
                      phi = array(1, dim = c(4000, I, J))))
#>   budget lambda1 lambda2 K       N   Utility
#> 1  1e+05    0.01    5000 1 4500000  9.990594
#> 2  1e+05    0.01    5000 2 2000000 14.971293
#> 3  1e+05    0.01    5000 3 1166667 17.508238
#> 4  1e+05    0.01    5000 4  750000 18.746519
#> 5  1e+05    0.01    5000 5  500000 19.383884
# }
```
