# occumb: Site Occupancy Modeling For Environmental DNA Metabarcoding

<!-- badges: start -->
[![R-CMD-check](https://github.com/fukayak/occumb/workflows/R-CMD-check/badge.svg)](https://github.com/fukayak/occumb/actions)
<!-- badges: end -->

occumb is an R package that aims to apply [the multispecies site occupancy modeling for environmental DNA (eDNA) metabarcoding](https://doi.org/10.1111/2041-210X.13732) easily.

This package allows the users to fit the model with a fully Bayesian approach, using the conventional formulas in R. It enables the analysis of the detectability of species in different stages of the workflow of eDNA metabarcoding and inference of species site occupancy while accounting for false negatives. It also provides the functionality for a model-based inference to assist the optimization of the study design.

See the [original paper](https://doi.org/10.1111/2041-210X.13732) for details of the model and inference.

## Important note
The current programs are an alpha version released to present the concept of package development. In other words, please be aware of the following points when you want to try this package:

- Only a few planned features have been implemented.
- The functions may not have been thoroughly tested.
- Documentation may be insufficient.
- Specifications may change significantly in the future.

## Installation

You can install the released version of occumb from the GitHub repository using the `install_github()` function in remote (or devtools) package:

``` r
remote::install_github("fukayak/occumb", ref = "main")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(occumb)
## basic example code
```

