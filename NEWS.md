# occumb (development version)
* `plot` and `summary` methods are exported correctly.

# occumb 1.0.0 (2023/09/14)
* Add `z`, `theta`, and `phi` arguments to `eval_util_L()`.
* Add `psi`, `theta`, and `phi` arguments to `eval_util_R()`.
* Fix `eval_util_R()` to account for site-heterogeneity of parameters.
* Remove `loglik()` function from the package.
* Some bug fixes, internal changes, and document improvements.

# occumb 0.6.1 (2023/07/27)
* Add model specification vignette.
* Some fixes and improvements of documentation.

# occumb v0.6.0 (2023/07/20)
* Add package vignette.
* Add `fish` and `fish_raw` data.
* Add `occumbGof` class.
* Add and fix methods for `occumbData`, `occumbFit`, and `occumbGof` classes.
* Change defaults for `cores` arguments for `gof()`, `eval_util_L()`, and `eval_util_R()` functions.
* Some bug fixes, internal changes, and documentation improvements.

# occumb v0.5.1 (2023/04/27)
* Add pkgdown website.
* Some bug fixes and internal changes.

# occumb v0.5.0 (2023/04/26)
* Add `get_post_samples()` and `get_post_summary()` functions.
* Add option for parallel computation in `gof()` function.
* Fix some bugs in `occumb()` functions.
* Fix an issue of parallel computation of `eval_util_L()` and `eval_util_R()` functions on Windows.

# occumb v0.4.2 (2022/12/13)
* Fix an issue of parallel computation of `eval_util_L()` and `eval_util_R()` functions on Windows.

# occumb v0.4.1 (2022/11/16)
* Fix some computational issues in `eval_util_L()` and `eval_util_R()` functions.

# occumb v0.4.0 (2022/06/03)
* Add `eval_util_L()` and `eval_util_R()` functions.
* Add `list_cond_L()` and `list_cond_R()` functions.

# occumb v0.3.0 (2022/03/31)
* Add `gof()` function.
* Add `loglik()` function.
* Add `...` argument to `occumb()` function.
* Add `data` field in `occumbFit` class.

# occumb v0.2.1 (2021/11/21)
* A few bug fixes.
* Add validations for the inputs of `occumb()`.

# occumb v0.2.0 (2021/11/19)
* Change occumbData class specification and `occumbData()` function.
* Add `occumb()` function.
* Add methods for `occumbFit` class: `plot`, `print`, `summary`.

# occumb v0.1.0 (2021/9/10)
* Initial development of occumb package.
* Add `occumbData()` function.

