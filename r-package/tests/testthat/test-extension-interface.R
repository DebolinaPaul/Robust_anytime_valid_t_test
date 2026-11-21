test_that('test paths retain running validity ordering and stopping indices', {
  x <- c(-.3, .8, 1.2, -.4, 1.1, .6)
  prior <- ev_prior('gaussian', location = 3, scale = 1, alternative = 'greater')
  fit <- ev_ttest(x, prior, alpha = .2, tolerance = 2e-4)
  expect_s3_class(fit, 'ev_ttest')
  expect_equal(fit$path$n, 2:length(x))
  expect_true(all(diff(fit$path$p_eville) <= 2e-8))
  expect_true(all(diff(fit$path$p_ville) <= 1e-12))
  expect_true(all(fit$path$p_eville <= fit$path$p_ville + 1e-12))
  expect_equal(as.data.frame(fit), fit$path)
  for (method in c('eVille', 'Ville')) {
    values <- fit$path[[if (method == 'eVille') 'p_eville' else 'p_ville']]
    hits <- which(values <= .2)
    expected <- if (length(hits)) fit$path$n[hits[1]] else NA_integer_
    expect_equal(unname(fit$stopping_times[method]), expected)
  }
  wrong <- ev_calibration(ev_prior('gaussian', location = 3, scale = 1,
                                    alternative = 'less'))
  expect_error(ev_ttest(x, prior, calibration = wrong), 'match')
  expect_error(ev_ttest(x, prior, n0 = 3), 'n0')
  expect_error(ev_ttest(x, prior, alpha = 1), 'alpha')
  expect_error(ev_ttest(c(1, 1, 2), prior), 'variance')
})

test_that('one-sided Lebesgue paths have no Ville comparison', {
  x <- c(-.3, .8, 1.2, -.4, 1.1, .6)
  p <- ev_prior('lebesgue', alternative = 'greater')
  fit <- ev_ttest(x, p)
  expect_true(all(is.na(fit$path$p_ville)))
  expect_true(is.na(fit$stopping_times['Ville']))
  expect_true(all(diff(fit$path$p_eville) <= 1e-12))
  expect_true(all(fit$path$calibration_gap == 0))
})
