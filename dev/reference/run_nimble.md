# Run MCMC with the NIMBLE backend (internal)

Internal worker used to fit the model with the NIMBLE engine. This
function builds a NIMBLE model, runs MCMC (optionally in parallel),
summarizes the posterior samples, and returns an object compatible with
the jagsUI-style fit object.

## Usage

``` r
run_nimble(
  data,
  const,
  inits,
  params,
  model_code_strings,
  model_file,
  n.chains,
  n.iter,
  n.burnin,
  n.thin,
  parallel,
  ...
)
```

## Arguments

- data:

  A named list containing the observed data and covariates for NIMBLE.

- const:

  A named list of model constants for NIMBLE.

- inits:

  A function that returns a named list of initial values for the model
  parameters (one set of initial values will be generated per chain).

- params:

  Character vector of node names to monitor.

- model_code_strings:

  Character vector of NIMBLE model code lines to be converted to
  `nimbleCode`.

- model_file:

  A string identifying the model file (kept for compatibility with other
  engines); stored in the returned object.

- n.chains:

  Number of MCMC chains.

- n.iter:

  Total number of MCMC iterations per chain.

- n.burnin:

  Number of burn-in iterations.

- n.thin:

  Thinning interval.

- parallel:

  Logical; if `TRUE`, run chains in parallel using the parallel package.

- ...:

  Additional control arguments:

  `seed`

  :   See `nimble::runMCMC(setSeed = ...)`. `FALSE`: no seeding; `TRUE`:
      seed chain i with i; numeric vector (`length = n.chains`):
      per-chain seeds.

  `n.cores`

  :   Number of worker processes when `parallel=TRUE`. Defaults to
      [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html).

  `store.data`

  :   Logical; if `TRUE`, store the input `data` and generated initial
      values in the returned object.

  `verbose`

  :   Logical; if not `NULL`, temporarily set `nimble` options `verbose`
      and `MCMCprogressBar` accordingly during execution.

## Value

jagsUI-style fit object.
