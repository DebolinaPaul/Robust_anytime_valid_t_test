test_that("two-sided Rouder evidence agrees with BayesFactor", {
  skip_if_not_installed("BayesFactor")
  x <- c(-0.4, 0.9, 0.2, 1.1, -0.7, 0.5, 0.3, 1.4,
         -0.2, 0.8, 0.1, 0.6)
  for (r in c(1 / sqrt(2), 1)) {
    for (n in c(2, 5, 12)) {
      xx <- x[seq_len(n)]
      bf <- BayesFactor::ttestBF(x = xx, mu = 0, rscale = r,
                                nullInterval = NULL, posterior = FALSE)
      log_bf <- as.numeric(BayesFactor::extractBF(
        bf, logbf = TRUE, onlybf = TRUE))
      expect_lt(abs(ev_log_evidence(xx, "rouder", r) - log_bf), 2e-6)
    }
  }
})
