# Find sampler indices (fast)

This function provides a fast alternative to
`conf$findSamplersOnNodes()`, which can be slow for large target
samplers.

## Usage

``` r
find_sampler_indices_fast(conf, nodes)
```

## Arguments

- conf:

  A [`MCMCconf`](https://rdrr.io/pkg/nimble/man/MCMCconf-class.html)
  object from the nimble package.

- nodes:

  A parameter name (character string).

## Value

An integer vector of indices into `conf$getSamplers()`.
