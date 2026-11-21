test_that('a bracket crossing alpha is reported rather than silently resolved', {
  x <- c(-.3, .8)
  p <- ev_prior('gaussian', location = 3, scale = .7, alternative = 'greater')
  cal <- ev_calibration(p, method = "deterministic", tolerance = 1e-3)
  bracket <- ev_pvalues(ev_log_mixture(x, p), cal)
  expect_gt(bracket$gap, 0)
  alpha <- (bracket$lower + bracket$p_eville) / 2
  expect_warning(fit <- ev_ttest(x, p, alpha = alpha, calibration = cal),
                 'bracket straddles alpha')
  expect_equal(fit$calibration_ambiguous_n, 2L)
  expect_true(is.na(fit$stopping_times['eVille']))
})
