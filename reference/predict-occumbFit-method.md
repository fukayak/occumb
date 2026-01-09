# Predict method for occumbFit class.

Obtain predictions of parameters related to species occupancy and
detection from an `occumbFit` model object.

## Usage

``` r
# S4 method for class 'occumbFit'
predict(
  object,
  newdata = NULL,
  parameter = c("phi", "theta", "psi"),
  scale = c("response", "link"),
  type = c("quantiles", "mean", "samples"),
  output_dataframe = FALSE
)
```

## Arguments

- object:

  An `occumbFit` object.

- newdata:

  An optional `occumbData` object with covariates to be used for
  prediction. If omitted, the fitted covariates are used.

- parameter:

  The parameter to be predicted.

- scale:

  The scale on which the prediction is made. `type = "response"` returns
  the prediction on the original scale of the parameter. `type = "link"`
  returns the prediction on the link scale of the parameter.

- type:

  The type of prediction. `type = "quantiles"` returns 50% quantile as
  the posterior median of the prediction in addition to 2.5 and 97.5%
  quantiles as the lower and upper limits of the 95% credible interval
  of the prediction. `type = "mean"` returns the posterior mean of the
  prediction. `type = "samples"` returns the posterior samples of the
  prediction.

- output_dataframe:

  If `TRUE`, results are returned in data frame format.

## Value

Predictions are obtained as a matrix or array that can have dimensions
corresponding to statistics (or samples), species, sites, and
replicates. The `dimension` and `label` attributes are added to the
output object to inform these dimensions. If the sequence read count
data `y` has species, site, or replicate names appended as the
`dimnames` attribute (see Details in
[`occumbData()`](https://fukayak.github.io/occumb/reference/occumbData.md)),
they will be copied into the `label` attribute of the returned object.

When `output_dataframe = TRUE`, the results are returned in data frame
format where the attributes obtained when `output_dataframe = FALSE` are
incorporated into the table.

## Details

Applying [`predict()`](https://rdrr.io/r/stats/predict.html) to an
`occumbFit` object generates predictions for the specified parameter
(`phi`, `theta`, or `psi`) based on the estimated effects and the given
covariates. It is important to recognize that the predictions are
specific to the individual species being modeled since they depend on
the estimated species-specific effects (i.e., `alpha`, `beta`, and
`gamma`; see [the package
vignette](https://fukayak.github.io/occumb/articles/model_specification.html)
for details). When providing `newdata`, it must thus be assumed that the
set of species contained in `newdata` is the same as that of the data
being fitted.
