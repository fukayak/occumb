# Plot method for occumbFit class.

Applies [jagsUI](https://cran.r-project.org/package=jagsUI)'s plot
method to an `occumbFit` object to draw trace plots and density plots of
MCMC samples of model parameters.

## Usage

``` r
# S4 method for class 'occumbFit'
plot(x, y = NULL, ...)
```

## Arguments

- x:

  An `occumbFit` object.

- y:

  `NULL`

- ...:

  Additional arguments passed to the plot method for
  [jagsUI](https://cran.r-project.org/package=jagsUI) object.

## Value

Returns `NULL` invisibly.
