test_that("Lebesgue calibration matches analytic values", {
  expect_equal(
    g_lebesgue(c(0, sqrt(pi), Inf)),
    c(1, 1, 0)
  )

  expect_equal(
    g_lebesgue(2 * sqrt(pi)),
    0.5 + 1 / pi,
    tolerance = 1e-12
  )
})

