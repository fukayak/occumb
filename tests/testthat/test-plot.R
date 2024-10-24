### Test for plot() occumbFit --------------------------------------------------
test_that("plot() has known output for occumbFit)", {
  plot_occumbFit <- function() plot(occumb:::internal_fit)
  vdiffr::expect_doppelganger(
    title = "occumbFit",
    fig = plot_occumbFit
  )
})

### Test for plot() gof --------------------------------------------------------
test_that("plot() has known output for gof", {
  plot_gof <- function() plot(occumb:::gof_ft)
  vdiffr::expect_doppelganger(
    title = "gof",
    fig = plot_gof
  )
})
