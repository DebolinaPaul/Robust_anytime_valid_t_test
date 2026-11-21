# Independent integral: integrate over the chi_n radial variable, rather
# than evaluate the closed-form mixture used by the production functions.
radial_reference <- function(x, prior, scale = 1) {
  n <- length(x)
  u2 <- sum(x)^2 / sum(x^2)
  if (prior == "gaussian") {
    normalizer <- 1 + n * scale^2
    coefficient <- scale^2 * u2 / normalizer
    log_factor <- -0.5 * log(normalizer)
  } else {
    coefficient <- u2 / n
    log_factor <- 0.5 * log(2 * pi / n)
  }
  integrand <- function(r) {
    out <- numeric(length(r))
    ok <- r > 0
    out[ok] <- exp(log(2 * r[ok]) +
      stats::dchisq(r[ok]^2, df = n, log = TRUE) +
      coefficient * r[ok]^2 / 2 + log_factor)
    out
  }
  log(stats::integrate(integrand, 0, Inf, rel.tol = 1e-10,
                       abs.tol = 1e-11, subdivisions = 500L)$value)
}

test_that("Gaussian and Lebesgue formulas agree with radial integration", {
  cases <- list(c(-0.4, 0.9), c(-1, 0, 0.5, 1.2, -0.3),
                c(0.5, 1.2, -0.1, 0.8, 0.3, 1.1))
  for (x in cases) {
    for (s in c(0.2, 1, 3)) {
      expect_lt(abs(ev_log_evidence(x, "gaussian", s) -
                    radial_reference(x, "gaussian", s)), 2e-9)
    }
    expect_lt(abs(ev_log_evidence(x, "lebesgue") -
                  radial_reference(x, "lebesgue")), 2e-9)
  }
})

test_that("evidence agrees with the existing Python engine", {
  x <- utils::read.csv(test_path("fixtures", "python_observations.csv"))$x
  ref <- utils::read.csv(test_path("fixtures", "python_evidence.csv"))
  errors <- vapply(seq_len(nrow(ref)), function(i) {
    row <- ref[i, ]
    xx <- if (row$case == "prefix") x[seq_len(row$n)] else
      seq_len(row$n) - (row$n + 1) / 2
    prior <- if (row$prior == "cauchy") "rouder" else row$prior
    scale <- if (prior == "lebesgue") NULL else row$scale
    ev_log_evidence(xx, prior, scale) - row$log_m
  }, numeric(1))
  expect_lt(max(abs(errors)), 2e-8)
})

test_that("positive rescaling and sign reversal preserve symmetric evidence", {
  x <- c(-0.4, 0.9, 0.2, 1.1, -0.7, 0.5, 0.3, 1.4)
  for (prior in c("gaussian", "rouder", "lebesgue")) {
    reference <- ev_log_evidence(x, prior)
    for (a in c(1e-100, 0.1, 10, 1e100)) {
      expect_lt(abs(ev_log_evidence(a * x, prior) - reference), 2e-9)
    }
    expect_lt(abs(ev_log_evidence(-x, prior) - reference), 2e-9)
    # For nonzero null mean, rescale/translate that mean as well.
    expect_lt(abs(ev_log_evidence(7 + 3 * x, prior, mu = 7) - reference), 2e-9)
  }
  original <- ev_path(x, "gaussian", scale = 30)
  rescaled <- ev_path(100 * x, "gaussian", scale = 30)
  expect_equal(rescaled, original, tolerance = 1e-10)
})

test_that("undefined observations fail explicitly", {
  expect_error(ev_log_evidence(c(1, NA)), "finite")
  expect_error(ev_log_evidence(c(1, Inf)), "finite")
  expect_error(ev_log_evidence(1), "length")
  expect_error(ev_log_evidence(c(1, 1)), "variance")
  expect_error(ev_log_evidence(c(0, 0)), "variance")
  expect_error(ev_log_evidence(c(1, 2), scale = -1), "scale")
  expect_error(ev_path(c(1, 1, 2)), "variance")
})
