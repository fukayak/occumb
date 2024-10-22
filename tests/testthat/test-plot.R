### Test for plot() occumbFit --------------------------------------------------
test_that("plot() works for occumbFit as expected", {
  vdiffr::expect_doppelganger(
    title = "occumbFit",
    fig = plot(internal_fit)
  )
})

### Test for plot() gof --------------------------------------------------------
test_that("plot() works for gof as expected", {
  vdiffr::expect_doppelganger(
    title = "gof",
    fig = plot(gof_ft)
  )
})
