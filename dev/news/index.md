# Changelog

## occumb (development version)

- Fix timezone specification for snapshot testing
- Add an option to select the NIMBLE engine in
  [`occumb()`](https://fukayak.github.io/occumb/dev/reference/occumb.md)

### occumb 1.2.2 (2026/1/9)

- This patch release includes internal changes for the anticipated JAGS
  5.0.0 release and fixes for math rendering issues in the package
  vignettes.

### occumb 1.2.1 (2025/7/16)

- This patch release only fixes bibliographic information in package
  documents.

### occumb 1.2.0 (2025/5/23)

- Add a data frame input option in
  [`occumbData()`](https://fukayak.github.io/occumb/dev/reference/occumbData.md).
- Add a data frame output option in
  [`get_post_samples()`](https://fukayak.github.io/occumb/dev/reference/get_posterior.md)
  and
  [`get_post_summary()`](https://fukayak.github.io/occumb/dev/reference/get_posterior.md).
- Add a data frame output option in
  [`predict()`](https://rdrr.io/r/stats/predict.html).
- Add `stats = chi_squared` option to
  [`gof()`](https://fukayak.github.io/occumb/dev/reference/gof.md).
- Internal changes to fix a number of known bugs, more helpful error
  messages, and additional testing.
- Improved function documentation.

### occumb 1.1.0 (2024/3/26)

- Add [`predict()`](https://rdrr.io/r/stats/predict.html) method for
  `occumbFit` class.
- Fix [`summary()`](https://rdrr.io/r/base/summary.html) method for
  `occumbFit`: no longer outputs comments on convergence and DIC.
- Internal changes to fix a number of known bugs, more helpful messages,
  and additional testing.
- Improved function documentation and vignette.

### occumb 1.0.3 (2024/01/04)

- Internal changes to fix a number of known bugs, more helpful error
  messages, and additional testing.
- Improved function documentation.

### occumb 1.0.2 (2023/10/19)

- Fix license issue: occumb is licensed under GPLv3.
- Some document fixes.

### occumb 1.0.1 (2023/09/21)

- This patch release fixes issues with
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods not being
  exported correctly.
- It also improves
  [`gof()`](https://fukayak.github.io/occumb/dev/reference/gof.md) to
  accept additional arguments for figure formatting.

### occumb 1.0.0 (2023/09/14)

- Add `z`, `theta`, and `phi` arguments to
  [`eval_util_L()`](https://fukayak.github.io/occumb/dev/reference/eval_util_L.md).
- Add `psi`, `theta`, and `phi` arguments to
  [`eval_util_R()`](https://fukayak.github.io/occumb/dev/reference/eval_util_R.md).
- Fix
  [`eval_util_R()`](https://fukayak.github.io/occumb/dev/reference/eval_util_R.md)
  to account for site-heterogeneity of parameters.
- Remove `loglik()` function from the package.
- Some bug fixes, internal changes, and document improvements.

### occumb 0.6.1 (2023/07/27)

- Add model specification vignette.
- Some fixes and improvements of documentation.

### occumb 0.6.0 (2023/07/20)

- Add package vignette.
- Add `fish` and `fish_raw` data.
- Add `occumbGof` class.
- Add and fix methods for `occumbData`, `occumbFit`, and `occumbGof`
  classes.
- Change defaults for `cores` arguments for
  [`gof()`](https://fukayak.github.io/occumb/dev/reference/gof.md),
  [`eval_util_L()`](https://fukayak.github.io/occumb/dev/reference/eval_util_L.md),
  and
  [`eval_util_R()`](https://fukayak.github.io/occumb/dev/reference/eval_util_R.md)
  functions.
- Some bug fixes, internal changes, and documentation improvements.

### occumb 0.5.1 (2023/04/27)

- Add pkgdown website.
- Some bug fixes and internal changes.

### occumb 0.5.0 (2023/04/26)

- Add
  [`get_post_samples()`](https://fukayak.github.io/occumb/dev/reference/get_posterior.md)
  and
  [`get_post_summary()`](https://fukayak.github.io/occumb/dev/reference/get_posterior.md)
  functions.
- Add option for parallel computation in
  [`gof()`](https://fukayak.github.io/occumb/dev/reference/gof.md)
  function.
- Fix some bugs in
  [`occumb()`](https://fukayak.github.io/occumb/dev/reference/occumb.md)
  functions.
- Fix an issue of parallel computation of
  [`eval_util_L()`](https://fukayak.github.io/occumb/dev/reference/eval_util_L.md)
  and
  [`eval_util_R()`](https://fukayak.github.io/occumb/dev/reference/eval_util_R.md)
  functions on Windows.

### occumb 0.4.2 (2022/12/13)

- Fix an issue of parallel computation of
  [`eval_util_L()`](https://fukayak.github.io/occumb/dev/reference/eval_util_L.md)
  and
  [`eval_util_R()`](https://fukayak.github.io/occumb/dev/reference/eval_util_R.md)
  functions on Windows.

### occumb 0.4.1 (2022/11/16)

- Fix some computational issues in
  [`eval_util_L()`](https://fukayak.github.io/occumb/dev/reference/eval_util_L.md)
  and
  [`eval_util_R()`](https://fukayak.github.io/occumb/dev/reference/eval_util_R.md)
  functions.

### occumb 0.4.0 (2022/06/03)

- Add
  [`eval_util_L()`](https://fukayak.github.io/occumb/dev/reference/eval_util_L.md)
  and
  [`eval_util_R()`](https://fukayak.github.io/occumb/dev/reference/eval_util_R.md)
  functions.
- Add
  [`list_cond_L()`](https://fukayak.github.io/occumb/dev/reference/list_cond_L.md)
  and
  [`list_cond_R()`](https://fukayak.github.io/occumb/dev/reference/list_cond_R.md)
  functions.

### occumb 0.3.0 (2022/03/31)

- Add [`gof()`](https://fukayak.github.io/occumb/dev/reference/gof.md)
  function.
- Add `loglik()` function.
- Add `...` argument to
  [`occumb()`](https://fukayak.github.io/occumb/dev/reference/occumb.md)
  function.
- Add `data` field in `occumbFit` class.

### occumb 0.2.1 (2021/11/21)

- A few bug fixes.
- Add validations for the inputs of
  [`occumb()`](https://fukayak.github.io/occumb/dev/reference/occumb.md).

### occumb 0.2.0 (2021/11/19)

- Change occumbData class specification and
  [`occumbData()`](https://fukayak.github.io/occumb/dev/reference/occumbData.md)
  function.
- Add
  [`occumb()`](https://fukayak.github.io/occumb/dev/reference/occumb.md)
  function.
- Add methods for `occumbFit` class: `plot`, `print`, `summary`.

### occumb 0.1.0 (2021/9/10)

- Initial development of occumb package.
- Add
  [`occumbData()`](https://fukayak.github.io/occumb/dev/reference/occumbData.md)
  function.
