# Package index

## Setup dataset

Functions for building a dataset object used for model-fitting.

- [`occumbData()`](https://fukayak.github.io/occumb/dev/reference/occumbData.md)
  : Constructor for occumbData data class.

## Model fitting

Functions for fitting a model.

- [`occumb()`](https://fukayak.github.io/occumb/dev/reference/occumb.md)
  : Model-fitting function.

## Model assessment

Functions for assessment of the fitted model.

- [`get_post_samples()`](https://fukayak.github.io/occumb/dev/reference/get_posterior.md)
  [`get_post_summary()`](https://fukayak.github.io/occumb/dev/reference/get_posterior.md)
  : Extract posterior samples or summary of parameters from a model-fit
  object.
- [`gof()`](https://fukayak.github.io/occumb/dev/reference/gof.md) :
  Goodness-of-fit assessment of the fitted model.

## Study design for effective species detection

Functions for the model-based assessment of study design.

- [`eval_util_L()`](https://fukayak.github.io/occumb/dev/reference/eval_util_L.md)
  : Expected utility for local species diversity assessments.
- [`eval_util_R()`](https://fukayak.github.io/occumb/dev/reference/eval_util_R.md)
  : Expected utility for regional species diversity assessments.
- [`list_cond_L()`](https://fukayak.github.io/occumb/dev/reference/list_cond_L.md)
  : Conditions for local assessment under certain budget and cost
  values.
- [`list_cond_R()`](https://fukayak.github.io/occumb/dev/reference/list_cond_R.md)
  : Conditions for regional assessment under certain budget and cost
  values.

## S4 methods

Generic function methods applied to occumb package outputs.

- [`plot(`*`<occumbFit>`*`)`](https://fukayak.github.io/occumb/dev/reference/plot-occumbFit-method.md)
  : Plot method for occumbFit class.
- [`plot(`*`<occumbGof>`*`)`](https://fukayak.github.io/occumb/dev/reference/plot-occumbGof-method.md)
  : Plot method for occumbGof class.
- [`predict(`*`<occumbFit>`*`)`](https://fukayak.github.io/occumb/dev/reference/predict-occumbFit-method.md)
  : Predict method for occumbFit class.
- [`summary(`*`<occumbData>`*`)`](https://fukayak.github.io/occumb/dev/reference/summary-occumbData-method.md)
  : Summary method for occumbData class.
- [`summary(`*`<occumbFit>`*`)`](https://fukayak.github.io/occumb/dev/reference/summary-occumbFit-method.md)
  : Summary method for occumbFit class.

## Dataset

Example dataset to illustrate function usage.

- [`fish`](https://fukayak.github.io/occumb/dev/reference/fish.md) :
  Fish eDNA metabarcoding dataset
- [`fish_raw`](https://fukayak.github.io/occumb/dev/reference/fish_raw.md)
  : Fish eDNA metabarcoding dataset
