# Extract posterior samples or summary of parameters from a model-fit object.

`get_post_samples()` extracts posterior samples of the specified
parameters from a model-fit object.

`get_post_summary()` extracts posterior summary of the specified
parameters from a model-fit object.

## Usage

``` r
get_post_samples(
  fit,
  parameter = c("z", "pi", "phi", "theta", "psi", "alpha", "beta", "gamma",
    "alpha_shared", "beta_shared", "gamma_shared", "Mu", "sigma", "rho"),
  output_dataframe = FALSE
)

get_post_summary(
  fit,
  parameter = c("z", "pi", "phi", "theta", "psi", "alpha", "beta", "gamma",
    "alpha_shared", "beta_shared", "gamma_shared", "Mu", "sigma", "rho"),
  output_dataframe = FALSE
)
```

## Arguments

- fit:

  An `occumbFit` object.

- parameter:

  A string of parameter name. See Details for possible choices and
  corresponding parameters.

- output_dataframe:

  If `TRUE`, results are returned in data frame format.

## Value

By default, `get_post_samples()` returns a vector, matrix, or array of
posterior samples for a selected parameter.

`get_post_summary()` returns, by default, a table (matrix) of the
posterior summary of the selected parameters. The elements of the
posterior summary are the same as those obtained with the
[`jags()`](https://kenkellner.com/jagsUI/reference/jags.html) function
in the `jagsUI` package: they include the mean, standard deviation,
percentiles of posterior samples; the `Rhat` statistic; the effective
sample size, `n.eff`; `overlap0`, which checks if 0 falls in the
parameter's 95% credible interval; and the proportion of the posterior
with the same sign as the mean, `f`.

The `dimension` and `label` attributes of the output object provide
information regarding the dimensions of the parameter.

When `output_dataframe = TRUE`, the results are returned in data frame
format where the attributes obtained when `output_dataframe = FALSE` are
incorporated into the table.

## Details

The functions return posterior samples or a summary of one of the
following parameters in the model, stored in the model-fit object `fit`:

- `z`:

  Site occupancy status of species.

- `pi`:

  Multinomial probabilities of species sequence read counts.

- `phi`:

  Sequence relative dominance of species.

- `theta`:

  Sequence capture probabilities of species.

- `psi`:

  Site occupancy probabilities of species.

- `alpha`:

  Species-specific effects on sequence relative dominance (`phi`).

- `beta`:

  Species-specific effects on sequence capture probabilities (`theta`).

- `gamma`:

  Species-specific effects on site occupancy probabilities (`psi`).

- `alpha_shared`:

  Effects on sequence relative dominance (`phi`) common across species.

- `beta_shared`:

  Effects on sequence capture probabilities (`theta`) that are common
  across species.

- `gamma_shared`:

  Effects on site occupancy probabilities (`psi`) that are common across
  species.

- `Mu`:

  Community-level averages of species-specific effects (`alpha`, `beta`,
  `gamma`).

- `sigma`:

  Standard deviations of species-specific effects (`alpha`, `beta`,
  `gamma`).

- `rho`:

  Correlation coefficients of the species-specific effects (`alpha`,
  `beta`, `gamma`).

See [the package
vignette](https://fukayak.github.io/occumb/articles/model_specification.html)
for details of these parameters.

The parameter may have dimensions corresponding to species, sites,
replicates, and effects (covariates), and when
`output_dataframe = FALSE`, the `dimension` and `label` attributes are
added to the output object to inform these dimensions. If the sequence
read count data `y` have species, site, or replicate names appended as
the `dimnames` attribute (see Details in
[`occumbData()`](https://fukayak.github.io/occumb/reference/occumbData.md)),
they are copied into the `label` attribute of the returned object.

## Examples

``` r
# \donttest{
# Generate the smallest random dataset (2 species * 2 sites * 2 reps)
I <- 2 # Number of species
J <- 2 # Number of sites
K <- 2 # Number of replicates
y_named <- array(sample.int(I * J * K), dim = c(I, J, K))
dimnames(y_named) <- list(c("species 1", "species 2"),
                          c("site 1", "site 2"), NULL)
data_named <- occumbData(y = y_named)

# Fitting a null model
fit <- occumb(data = data_named, n.iter = 10100)
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
#> Sampling from joint posterior, 100 iterations x 4 chains 
#>  
#> 
#> Calculating statistics....... 
#> 
#> Done. 

# Extract posterior samples
(post_sample_z <- get_post_samples(fit, "z"))
#> , , 1
#> 
#>       [,1] [,2]
#>  [1,]    1    1
#>  [2,]    1    1
#>  [3,]    1    1
#>  [4,]    1    1
#>  [5,]    1    1
#>  [6,]    1    1
#>  [7,]    1    1
#>  [8,]    1    1
#>  [9,]    1    1
#> [10,]    1    1
#> [11,]    1    1
#> [12,]    1    1
#> [13,]    1    1
#> [14,]    1    1
#> [15,]    1    1
#> [16,]    1    1
#> [17,]    1    1
#> [18,]    1    1
#> [19,]    1    1
#> [20,]    1    1
#> [21,]    1    1
#> [22,]    1    1
#> [23,]    1    1
#> [24,]    1    1
#> [25,]    1    1
#> [26,]    1    1
#> [27,]    1    1
#> [28,]    1    1
#> [29,]    1    1
#> [30,]    1    1
#> [31,]    1    1
#> [32,]    1    1
#> [33,]    1    1
#> [34,]    1    1
#> [35,]    1    1
#> [36,]    1    1
#> [37,]    1    1
#> [38,]    1    1
#> [39,]    1    1
#> [40,]    1    1
#> 
#> , , 2
#> 
#>       [,1] [,2]
#>  [1,]    1    1
#>  [2,]    1    1
#>  [3,]    1    1
#>  [4,]    1    1
#>  [5,]    1    1
#>  [6,]    1    1
#>  [7,]    1    1
#>  [8,]    1    1
#>  [9,]    1    1
#> [10,]    1    1
#> [11,]    1    1
#> [12,]    1    1
#> [13,]    1    1
#> [14,]    1    1
#> [15,]    1    1
#> [16,]    1    1
#> [17,]    1    1
#> [18,]    1    1
#> [19,]    1    1
#> [20,]    1    1
#> [21,]    1    1
#> [22,]    1    1
#> [23,]    1    1
#> [24,]    1    1
#> [25,]    1    1
#> [26,]    1    1
#> [27,]    1    1
#> [28,]    1    1
#> [29,]    1    1
#> [30,]    1    1
#> [31,]    1    1
#> [32,]    1    1
#> [33,]    1    1
#> [34,]    1    1
#> [35,]    1    1
#> [36,]    1    1
#> [37,]    1    1
#> [38,]    1    1
#> [39,]    1    1
#> [40,]    1    1
#> 
#> attr(,"parameter")
#> [1] "z"
#> attr(,"dimension")
#> [1] "Sample"  "Species" "Site"   
#> attr(,"label")
#> attr(,"label")$Sample
#> NULL
#> 
#> attr(,"label")$Species
#> [1] "species 1" "species 2"
#> 
#> attr(,"label")$Site
#> [1] "site 1" "site 2"
#> 
# Look dimensions of the parameter
attributes(post_sample_z)
#> $dim
#> [1] 40  2  2
#> 
#> $parameter
#> [1] "z"
#> 
#> $dimension
#> [1] "Sample"  "Species" "Site"   
#> 
#> $label
#> $label$Sample
#> NULL
#> 
#> $label$Species
#> [1] "species 1" "species 2"
#> 
#> $label$Site
#> [1] "site 1" "site 2"
#> 
#> 

# Extract posterior summary
(post_summary_z <- get_post_summary(fit, "z"))
#>        mean sd 2.5% 25% 50% 75% 97.5% Rhat n.eff overlap0 f
#> z[1,1]    1  0    1   1   1   1     1   NA     1        0 1
#> z[2,1]    1  0    1   1   1   1     1   NA     1        0 1
#> z[1,2]    1  0    1   1   1   1     1   NA     1        0 1
#> z[2,2]    1  0    1   1   1   1     1   NA     1        0 1
#> attr(,"parameter")
#> [1] "z"
#> attr(,"dimension")
#> [1] "Species" "Site"   
#> attr(,"label")
#> attr(,"label")$Species
#> [1] "species 1" "species 2"
#> 
#> attr(,"label")$Site
#> [1] "site 1" "site 2"
#> 
# Look dimensions of the parameter
attributes(post_summary_z)
#> $dim
#> [1]  4 11
#> 
#> $dimnames
#> $dimnames[[1]]
#> [1] "z[1,1]" "z[2,1]" "z[1,2]" "z[2,2]"
#> 
#> $dimnames[[2]]
#>  [1] "mean"     "sd"       "2.5%"     "25%"      "50%"      "75%"     
#>  [7] "97.5%"    "Rhat"     "n.eff"    "overlap0" "f"       
#> 
#> 
#> $parameter
#> [1] "z"
#> 
#> $dimension
#> [1] "Species" "Site"   
#> 
#> $label
#> $label$Species
#> [1] "species 1" "species 2"
#> 
#> $label$Site
#> [1] "site 1" "site 2"
#> 
#> 
# }
```
