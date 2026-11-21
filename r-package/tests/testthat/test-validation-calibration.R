test_that("analytic calibrations agree with null-angle integration", {
  for (s in c(0.2, 1, 3, 30)) {
    a <- sqrt(1 + 2 * s^2)
    for (c in c(0.1, 0.5, 1, 2, 10, 100)) {
      f <- function(v) pmin(1, a / (1 + 2 * s^2 * cos(v)^2) / c)
      # Explicit split at the kink makes this an independent quadrature
      # comparison without asking an integrator to find a moving kink.
      cuts <- c(0, pi / 2)
      if (c > 1 / a && c < a) {
        root <- stats::uniroot(function(v) {
          a / (1 + 2 * s^2 * cos(v)^2) - c
        }, c(0, pi / 2), tol = 1e-12)$root
        cuts <- c(0, root, pi / 2)
      }
      reference <- (2 / pi) * sum(vapply(seq_len(length(cuts) - 1L),
        function(k) stats::integrate(f, cuts[k], cuts[k + 1L],
          rel.tol = 1e-10, abs.tol = 1e-11, subdivisions = 500L)$value,
        numeric(1)))
      expect_lt(abs(ev_calibrate(log(c), "gaussian", s) - reference), 2e-8)
    }
  }
  for (c in c(0.5, sqrt(pi), 2 * sqrt(pi), 10, 100)) {
    f <- function(v) pmin(1, sqrt(pi) / (c * cos(v)^2))
    reference <- (2 / pi) * stats::integrate(f, 0, pi / 2,
      rel.tol = 1e-10, abs.tol = 1e-11, subdivisions = 500L)$value
    expect_lt(abs(ev_calibrate(log(c), "lebesgue") - reference), 2e-8)
  }
  expect_equal(ev_calibrate(log(2 * sqrt(pi)), "lebesgue"),
               0.5 + 1 / pi, tolerance = 1e-12)
  # Large log evidence must not underflow prematurely to p = 0.
  small_p <- ev_calibrate(1000, "lebesgue")
  expect_gt(small_p, 0)
  expect_lt(abs(log(small_p) - (log(4 / pi) + log(pi) / 4 - 500)), 1e-12)
  # Rounding immediately above the Gaussian lower boundary.
  s <- 1.4273044337433236
  lc <- -0.5 * log1p(2 * s^2) + .Machine$double.eps / 2
  expect_true(is.finite(ev_calibrate(lc, "gaussian", s)))
})

test_that("eVille is no larger than Ville for the same proper prior", {
  lc <- seq(-10, 20, length.out = 301)
  for (s in c(0.2, 1, 3, 30)) {
    pe <- ev_calibrate(lc, "gaussian", s)
    pv <- ev_calibrate(lc, "gaussian", s, method = "Ville")
    expect_true(all(pe <= pv + 1e-12))
    expect_true(all(diff(pe) <= 1e-12))
    expect_true(all(pe >= 0 & pe <= 1 + 1e-12))
  }
  lc <- log(c(0.2, 1, 2, 5, 20))
  pe <- ev_calibrate(lc, "rouder")
  pv <- ev_calibrate(lc, "rouder", method = "Ville")
  expect_true(all(pe <= pv + 2e-8))
  expect_true(all(diff(pe) <= 2e-8))
  expect_true(all(pe >= 0 & pe <= 1 + 2e-8))
  # Independent clipped-angle quadrature at moderate thresholds.
  for (c in c(1, 2)) {
    f <- function(v) vapply(v, function(vv) {
      exp(min(0, .ev_log_rouder(2, 2 * cos(vv)^2, 1 / sqrt(2)) - log(c)))
    }, numeric(1))
    reference <- (2 / pi) * stats::integrate(f, 0, pi / 2,
      rel.tol = 1e-9, abs.tol = 0, subdivisions = 500L)$value
    expect_lt(abs(ev_calibrate(log(c), "rouder") - reference), 2e-8)
  }
  # This checks a strict improvement without requiring it on every dataset.
  expect_lt(ev_calibrate(log(20), "gaussian", 30), 0.05)
  expect_equal(ev_calibrate(log(20), "gaussian", 1), 0.05)
  # A very broad Cauchy prior has a narrow null-angle rejection region.
  # Direct unsplit angle integration can miss it entirely.
  broad <- ev_calibrate(0, "rouder", scale = 1e20)
  asymptotic <- 4 / (pi^(5 / 4) * sqrt(1e20))
  expect_lt(abs(broad / asymptotic - 1), 1e-7)
})

test_that("running paths start at two and omit Ville for Lebesgue", {
  x <- c(-0.4, 0.9, 0.2, 1.1, -0.7, 0.5)
  p <- ev_path(x, "gaussian", scale = 30)
  expect_equal(p$n, 2:length(x))
  expect_true(all(diff(p$p_eville) <= 1e-12))
  expect_true(all(diff(p$p_ville) <= 1e-12))
  expect_true(all(p$p_eville <= p$p_ville + 1e-12))
  flat <- ev_path(x, "lebesgue")
  expect_true(all(is.na(flat$p_ville)))
  expect_true(all(diff(flat$p_eville) <= 1e-12))
  expect_error(ev_calibrate(1, "lebesgue", method = "Ville"), "not defined")
  expect_equal(ev_calibrate(c(-Inf, Inf), "gaussian"), c(1, 0))
  expect_equal(ev_calibrate(c(-Inf, Inf), "lebesgue"), c(1, 0))
})
