# Conditions for local assessment under certain budget and cost values.

`list_cond_L()` constructs a list of possible local species diversity
assessment conditions under the specified budget and cost values.

## Usage

``` r
list_cond_L(budget, lambda1, lambda2, fit, K = NULL)
```

## Arguments

- budget:

  A numeric specifying budget amount. The currency unit is arbitrary but
  must be consistent with that of `lambda1` and `lambda2`.

- lambda1:

  A numeric specifying the cost per sequence read for high-throughput
  sequencing. The currency unit is arbitrary but must be consistent with
  that of `budget` and `lambda2`.

- lambda2:

  A numeric specifying the cost per replicate for library preparation.
  The currency unit is arbitrary but must be consistent with that of
  `budget` and `lambda1`.

- fit:

  An `occumbFit` object.

- K:

  An optional vector for manually specifying the number of replicates.

## Value

A data frame containing columns named `budget`, `lambda1`, `lambda2`,
`K`, and `N`.

## Details

This function can generate a data frame object to be given to the
`settings` argument of
[`eval_util_L()`](https://fukayak.github.io/occumb/reference/eval_util_L.md);
see Examples of
[`eval_util_L()`](https://fukayak.github.io/occumb/reference/eval_util_L.md).
By default, it outputs a list of all feasible combinations of values for
the number of replicates per site `K` and the sequencing depth per
replicate `N` based on the given budget, cost values, and number of
sites (identified by reference to the `fit` object). The resulting `N`
can be a non-integer because it is calculated simply by assuming that
the maximum value can be obtained. To obtain a list for only a subset of
the possible `K` values under a given budget and cost value, the `K`
argument is used to provide a vector of the desired `K` values.
