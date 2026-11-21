test_that('directional calibrations bracket the independent Python values', {
  ref <- utils::read.csv(test_path('fixtures', 'shifted_calibration.csv'))
  key <- paste(ref$prior, ref$location, ref$scale, ref$alternative)
  for (rows in split(seq_len(nrow(ref)), key)) {
    group <- ref[rows, ]
    prior <- extension_prior_row(group[1, ])
    cal <- ev_calibration(prior, method = "deterministic", tolerance = 1e-4)
    result <- ev_pvalues(group$log_c, cal)
    expect_true(all(result$lower <= group$g + 2e-7))
    expect_true(all(result$p_eville >= group$g - 2e-7))
    expect_true(all(result$gap <= 1e-4 + 1e-12))
    expect_equal(result$gap, result$p_eville - result$lower, tolerance = 1e-12)
    if (prior$proper) {
      expect_true(all(result$p_eville <= exp(pmin(0, -group$log_c)) + 1e-12))
    }
  }
})

test_that('convex cell bounds contain independently integrated clipped curves', {
  # Exact positive/negative exponential components provide a convex curve
  # with an interior minimum and an independently evaluated reference.
  at <- function(v) c(positive = 2 * sqrt(2) * sin(v),
                      negative = log(.3) - 3 * sqrt(2) * sin(v))
  for (ends in list(c(-1.3, -.4), c(-.6, .5), c(.2, 1.2))) {
    v <- c(ends[1], mean(ends), ends[2])
    cell <- list(v = v, parts = t(vapply(v, at, c(positive = 0, negative = 0))))
    for (threshold in c(.5, 2, 10)) {
      bounds <- .evx_cell_bounds(cell, log(threshold))
      fun <- function(z) pmin(1, (exp(2 * sqrt(2) * sin(z)) +
                                    .3 * exp(-3 * sqrt(2) * sin(z))) / threshold)
      reference <- stats::integrate(fun, ends[1], ends[2], rel.tol = 1e-10,
                                     subdivisions = 500L)$value / pi
      expect_lte(bounds['lower'], reference + 2e-9)
      expect_gte(bounds['upper'], reference - 2e-9)
    }
  }
})

test_that('shifted two-sided Gaussian calibration brackets a separate null integral', {
  p <- ev_prior('gaussian', location = 3, scale = .7)
  result <- ev_pvalues(0, ev_calibration(p, method = "deterministic", tolerance = 1e-4))
  # A shifted Gaussian mixture is a scaled noncentral-t density. This
  # independent oracle is used only at moderate noncentrality.
  D <- 1 + 2 * .7^2
  f <- function(v) {
    t <- tan(v)
    log_m <- stats::dt(t / sqrt(D), df = 1,
                      ncp = sqrt(2) * 3 / sqrt(D), log = TRUE) -
      .5 * log(D) - stats::dt(t, df = 1, log = TRUE)
    exp(pmin(0, log_m))
  }
  reference <- stats::integrate(f, -pi / 2, pi / 2, rel.tol = 1e-9,
                                subdivisions = 500L)$value / pi
  expect_lte(result$lower, reference + 2e-7)
  expect_gte(result$p_eville, reference - 2e-7)
  expect_lt(result$gap, 1e-4)
})

test_that('analytic calibrations preserve centered and full/half flat cases', {
  logs <- c(-Inf, -2, 0, 2, log(20), 1000, Inf)
  for (scale in c(1, 30)) {
    p <- ev_prior('gaussian', scale = scale)
    actual <- ev_pvalues(logs, ev_calibration(p, method = "deterministic"))
    expect_equal(actual$p_eville, ev_calibrate(logs, 'gaussian', scale),
                 tolerance = 1e-12)
    expect_true(all(actual$gap == 0))
  }
  full <- ev_pvalues(logs, ev_calibration(ev_prior('lebesgue')))$p_eville
  expect_equal(full, ev_calibrate(logs, 'lebesgue'), tolerance = 1e-12)
  greater <- ev_pvalues(logs, ev_calibration(ev_prior('lebesgue', alternative = 'greater')))
  less <- ev_pvalues(logs, ev_calibration(ev_prior('lebesgue', alternative = 'less')))
  expect_equal(greater, less)
  expect_lt(greater$p_eville[which(logs == 2)], full[which(logs == 2)])
  expect_gt(greater$p_eville[which(logs == 1000)], 0)
})
