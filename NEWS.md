# occumb (development version)
* Add test-coverage workflow
* Add snapshot tests for `plot()` and `summary()`
* Fix an issue with an error message by character covariates.
* Add `stats = chi_squared` option to `gof()`
* Fix an issue with having `.` with other covariates in formula
* Fix an issue with throwing a warning message instead of an error
* Change error message for validate_occumbData 
* Add a test to throw an error when integer is too large
* Add tests for `gof` under unbalanced designs
* Fix an issue with `sigma`'s default value.
* Add 'data.frame' as an input of occumbData.
* Fix a bug in attribute assignment in `get_post_samples` and `get_post_summary`.
* Fixed `get_post_samples()` not to output a vector.
* Add a data frame output option in `get_post_samples()` and `get_post_summary()`.

## occumb 1.1.0 (2024/3/26)
* Add `predict()` method for `occumbFit` class.
* Fix `summary()` method for `occumbFit`: no longer outputs comments on convergence and DIC.
* Internal changes to fix a number of known bugs, more helpful messages, and additional testing.
* Improved function documentation and vignette.

## occumb 1.0.3 (2024/01/04)
* Internal changes to fix a number of known bugs, more helpful error messages, and additional testing.
* Improved function documentation.

## occumb 1.0.2 (2023/10/19)
* Fix license issue: occumb is licensed under GPLv3.
* Some document fixes.

## occumb 1.0.1 (2023/09/21)
* This patch release fixes issues with `plot()` and `summary()` methods not being exported correctly.
* It also improves `gof()` to accept additional arguments for figure formatting.

## occumb 1.0.0 (2023/09/14)
* Add `z`, `theta`, and `phi` arguments to `eval_util_L()`.
* Add `psi`, `theta`, and `phi` arguments to `eval_util_R()`.
* Fix `eval_util_R()` to account for site-heterogeneity of parameters.
* Remove `loglik()` function from the package.
* Some bug fixes, internal changes, and document improvements.

## occumb 0.6.1 (2023/07/27)
* Add model specification vignette.
* Some fixes and improvements of documentation.

## occumb 0.6.0 (2023/07/20)
* Add package vignette.
* Add `fish` and `fish_raw` data.
* Add `occumbGof` class.
* Add and fix methods for `occumbData`, `occumbFit`, and `occumbGof` classes.
* Change defaults for `cores` arguments for `gof()`, `eval_util_L()`, and `eval_util_R()` functions.
* Some bug fixes, internal changes, and documentation improvements.

## occumb 0.5.1 (2023/04/27)
* Add pkgdown website.
* Some bug fixes and internal changes.

## occumb 0.5.0 (2023/04/26)
* Add `get_post_samples()` and `get_post_summary()` functions.
* Add option for parallel computation in `gof()` function.
* Fix some bugs in `occumb()` functions.
* Fix an issue of parallel computation of `eval_util_L()` and `eval_util_R()` functions on Windows.

## occumb 0.4.2 (2022/12/13)
* Fix an issue of parallel computation of `eval_util_L()` and `eval_util_R()` functions on Windows.

## occumb 0.4.1 (2022/11/16)
* Fix some computational issues in `eval_util_L()` and `eval_util_R()` functions.

## occumb 0.4.0 (2022/06/03)
* Add `eval_util_L()` and `eval_util_R()` functions.
* Add `list_cond_L()` and `list_cond_R()` functions.

## occumb 0.3.0 (2022/03/31)
* Add `gof()` function.
* Add `loglik()` function.
* Add `...` argument to `occumb()` function.
* Add `data` field in `occumbFit` class.

## occumb 0.2.1 (2021/11/21)
* A few bug fixes.
* Add validations for the inputs of `occumb()`.

## occumb 0.2.0 (2021/11/19)
* Change occumbData class specification and `occumbData()` function.
* Add `occumb()` function.
* Add methods for `occumbFit` class: `plot`, `print`, `summary`.

## occumb 0.1.0 (2021/9/10)
* Initial development of occumb package.
* Add `occumbData()` function.

