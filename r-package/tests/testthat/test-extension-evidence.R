test_that('all shifted and directional mixtures agree with Python references', {
  ref <- utils::read.csv(test_path('fixtures', 'shifted_evidence.csv'))
  error <- vapply(seq_len(nrow(ref)), function(i) {
    row <- ref[i, ]
    p <- extension_prior_row(row)
    ev_log_mixture(extension_data(row$n, row$t), p) - row$log_m
  }, numeric(1))
  expect_true(all(is.finite(error)))
  expect_lt(max(abs(error)), 1e-7)
})

test_that('one-sided priors are conditioned and reflection preserves evidence', {
  x <- extension_data(5, 1.7)
  for (family in c('gaussian', 'cauchy')) {
    full <- ev_prior(family, location = -2.3, scale = 1.1)
    positive <- ev_prior(family, location = -2.3, scale = 1.1,
                         alternative = 'greater')
    negative <- ev_prior(family, location = -2.3, scale = 1.1,
                         alternative = 'less')
    pos <- ev_log_mixture(x, positive)
    neg <- ev_log_mixture(x, negative)
    recombined <- .ev_logadd(full$log_positive_mass + pos,
                             full$log_negative_mass + neg)
    expect_lt(abs(recombined - ev_log_mixture(x, full)), 1e-8)
    reflected <- ev_prior(family, location = 2.3, scale = 1.1,
                          alternative = 'greater')
    expect_lt(abs(neg - ev_log_mixture(-x, reflected)), 1e-8)
    expect_lt(abs(ev_log_mixture(x, full) - ev_log_mixture(-x,
      ev_prior(family, location = 2.3, scale = 1.1))), 1e-8)
  }
  positive <- ev_log_mixture(x, ev_prior('lebesgue', alternative = 'greater'))
  negative <- ev_log_mixture(x, ev_prior('lebesgue', alternative = 'less'))
  expect_lt(abs(.ev_logadd(positive, negative) -
    ev_log_mixture(x, ev_prior('lebesgue'))), 1e-10)
})

test_that('rescaling preserves shifted and directional evidence', {
  x <- extension_data(5, -1.2)
  for (family in c('gaussian', 'cauchy', 'lebesgue')) {
    prior <- if (family == 'lebesgue') ev_prior(family, alternative = 'less') else
      ev_prior(family, location = -3, scale = .7, alternative = 'less')
    reference <- ev_log_mixture(x, prior)
    for (a in c(1e-80, .1, 20, 1e80)) {
      expect_lt(abs(ev_log_mixture(a * x, prior) - reference), 2e-8)
    }
    expect_lt(abs(ev_log_mixture(7 + 2 * x, prior, mu = 7) - reference), 2e-8)
  }
})

test_that('directional evidence obeys composite-null monotonicity', {
  x <- extension_data(6, -1)
  for (family in c('gaussian', 'cauchy', 'lebesgue')) {
    for (alternative in c('greater', 'less')) {
      p <- if (family == 'lebesgue') ev_prior(family, alternative = alternative) else
        ev_prior(family, location = 2, scale = 1, alternative = alternative)
      shifts <- seq(-1, 1, length.out = 5)
      value <- vapply(shifts, function(a) ev_log_mixture(x + a, p), numeric(1))
      if (alternative == 'greater') expect_true(all(diff(value) >= -1e-8))
      else expect_true(all(diff(value) <= 1e-8))
    }
  }
})

test_that('centered two-sided evidence retains the original reference cases', {
  x <- extension_data(8, 2.3)
  for (family in c('gaussian', 'cauchy', 'lebesgue')) {
    old <- if (family == 'cauchy') 'rouder' else family
    # The notebook's normal-scale quadrature uses relative tolerance 2e-9.
    # Its checked omitted tail can differ from the older centered-only
    # quadrature by about 1e-10; do not demand agreement beyond that accuracy.
    expect_equal(ev_log_mixture(x, ev_prior(family)),
                 ev_log_evidence(x, old), tolerance = 2e-9)
  }
  expect_error(ev_prior('rouder', location = 1), 'shifted Cauchy')
  expect_error(ev_prior('gaussian', scale = 0), 'scale')
  expect_error(ev_prior('lebesgue', scale = 1), 'Omit')
  edited <- ev_prior('gaussian')
  edited$location <- 1
  expect_error(ev_log_mixture(x, edited), 'modified')
})
