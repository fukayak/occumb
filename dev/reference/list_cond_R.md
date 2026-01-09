# Conditions for regional assessment under certain budget and cost values.

`list_cond_R()` constructs a list of possible regional species diversity
assessment conditions under the specified budget and cost values.

## Usage

``` r
list_cond_R(budget, lambda1, lambda2, lambda3, J = NULL, K = NULL)
```

## Arguments

- budget:

  A numeric specifying budget amount. The currency unit is arbitrary but
  must be consistent with that of `lambda1`, `lambda2`, and `lambda3`.

- lambda1:

  A numeric specifying the cost per sequence read for high-throughput
  sequencing. The currency unit is arbitrary but must be consistent with
  that of `budget`, `lambda2`, and `lambda3`.

- lambda2:

  A numeric specifying the cost per replicate for library preparation.
  The currency unit is arbitrary but must be consistent with that of
  `budget`, `lambda1`, and `lambda3`.

- lambda3:

  A numeric specifying the visiting cost per site. The currency unit is
  arbitrary but must be consistent with that of `budget`, `lambda1`, and
  `lambda2`.

- J:

  An optional vector for manually specifying the number of sites

- K:

  An optional vector used to specify the number of replicates manually.
  For computational convenience, the `K` values must be in ascending
  order.

## Value

A data frame containing columns named `budget`, `lambda1`, `lambda2`,
`lambda3`, `J`, `K`, and `N`.

## Details

This function can generate a data frame object to be given to the
`settings` argument of
[`eval_util_R()`](https://fukayak.github.io/occumb/dev/reference/eval_util_R.md);
see Examples of
[`eval_util_R()`](https://fukayak.github.io/occumb/dev/reference/eval_util_R.md).
By default, it outputs a list of all feasible combinations of values for
the number of sites `J`, number of replicates per site `K`, and
sequencing depth per replicate `N` based on the given budget and cost
values. The resulting `N` can be a non-integer because it is calculated
simply by assuming that the maximum value can be obtained. If one wants
to obtain a list for only a subset of the possible values of `J` and `K`
under a given budget and cost value, use the `J` and/or `K` arguments
(in fact, it is recommended that a relatively small number of `K` values
be specified using the `K` argument because the list of all conditions
achievable under moderate budget and cost values can be large, and it is
rarely practical to have a vast number of replicates per site). If a
given combination of `J` and `K` values is not feasible under the
specified budget and cost values, the combination will be ignored and
excluded from the output.
