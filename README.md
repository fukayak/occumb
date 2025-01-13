# occumb: Site Occupancy Modeling for Environmental DNA Metabarcoding

<!-- badges: start -->
[![R-CMD-check](https://github.com/fukayak/occumb/workflows/R-CMD-check/badge.svg)](https://github.com/fukayak/occumb/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/occumb)](https://CRAN.R-project.org/package=occumb)
[![Total downloads](http://cranlogs.r-pkg.org/badges/grand-total/occumb)](https://cran.r-project.org/package=occumb)
[![Monthly downloads](http://cranlogs.r-pkg.org/badges/last-month/occumb)](https://cran.r-project.org/package=occumb)
[![Codecov test coverage](https://codecov.io/gh/fukayak/occumb/graph/badge.svg)](https://app.codecov.io/gh/fukayak/occumb)
<!-- badges: end -->

occumb is an R package that provides functionalities for the application of [multispecies site occupancy modeling for environmental DNA (eDNA) metabarcoding](https://doi.org/10.1111/2041-210X.13732).

This package allows users to fit the model with a fully Bayesian approach using conventional formulas in R. This enables the analysis of the detectability of species at different stages of the workflow of eDNA metabarcoding and the inference of species site occupancy while accounting for false negatives. It also provides functionality for model-based inference to assist in the optimization of the study design.

See [the package vignette](https://fukayak.github.io/occumb/articles/occumb.html) to learn how to use the package and the [original paper](https://doi.org/10.1111/2041-210X.13732) for details of the model and inference.

## Installation
You need first install JAGS following the instructions on the [JAGS homepage](https://mcmc-jags.sourceforge.io/).

You can then install the latest stable version of the package from CRAN:

``` r
install.packages("occumb")
```

or the GitHub repository:

``` r
remotes::install_github("fukayak/occumb", ref = "main")
```

## Contact information

Questions and bug reports can be emailed to Keiichi Fukaya (fukaya.keiichi@nies.go.jp).

## Credits

The development of occumb would not have been possible without, among others, Martyn Plummer's [JAGS](https://mcmc-jags.sourceforge.io/) and Ken Kellner's [jagsUI](https://CRAN.R-project.org/package=jagsUI) R package, because the main functionality of occumb for model fitting via Markov chain Monte Carlo (MCMC) relies on these libraries. Taku Kadoya encouraged us to develop this package. Koji Makiyama, Shinya Uryu, and Kentaro Matsuura contributed to package testing. Mito Ikemoto and Kentaro Uehira provided valuable feedback through package trials and documentation reviews. Funding was provided by Japan Society for the Promotion of Science (KAKENHI; Nos. 20K06102 and 23H02240).

