# occumb: Site Occupancy Modeling For Environmental DNA Metabarcoding

<!-- badges: start -->
[![R-CMD-check](https://github.com/fukayak/occumb/workflows/R-CMD-check/badge.svg)](https://github.com/fukayak/occumb/actions)
<!-- badges: end -->

occumb is an R package that aims to apply [the multispecies site occupancy modeling for environmental DNA (eDNA) metabarcoding](https://doi.org/10.1111/2041-210X.13732) easily.

This package allows the users to fit the model, possibly with covariates, using a fully Bayesian approach to analyze the detectability of species in different stages of the workflow of eDNA metabarcoding and to infer species site occupancy while accounting for false negatives. It also provides the functionality for a model-based inference to assist the optimization of the study design.

See the [original paper](https://doi.org/10.1111/2041-210X.13732) for details of the model and inference.

## Important note
The current programs are an alpha version released to present the concept of package development. In other words, please be aware of the following points when you want to try this package:

- Only a limited number of planned features have been implemented.
- Functions may have not yet been thoroughly tested.
- Documentation may be insufficient.
- Specifications may be subject to significant changes in the future.

## Installation

You can install the released version of occumb from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("occumb")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(occumb)
## basic example code
```

