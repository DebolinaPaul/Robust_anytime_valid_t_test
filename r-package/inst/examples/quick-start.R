# One-sample t-test with the same Monte Carlo calibration as the notebook.
# Run after installing eVille. Preserve the stated observation order.
library(eVille)

x <- c(-0.3, -1.1, -0.7, -0.5, -1.4, -0.8, -0.2, -1.0)

# N(-4,1) conditioned on theta < 0; fixed before observing x.
prior <- ev_prior("gaussian", location = -4, scale = 1,
                  alternative = "less")
null_t <- ev_null_sample()
cal <- ev_calibration(prior, method = "monte_carlo", null_t = null_t)
fit <- ev_ttest(x, prior = prior, mu = 0, alpha = 0.05, calibration = cal)
print(summary(fit))
plot(fit, main = "One-sample t-test: Monte Carlo calibration")

# The numerical bracket does not account for Monte Carlo sampling error.
# Inspect empirical brackets and MC standard errors in the trajectory.
trajectory <- as.data.frame(fit)
print(tail(trajectory))
summary_table <- summary(fit)$results

# The half-Lebesgue benchmark has its own analytic calibration.
flat <- ev_ttest(x, prior = ev_prior("lebesgue", alternative = "less"),
                 mu = 0, alpha = 0.05)
print(summary(flat))

# Optional exports:
# write.csv(summary_table, "eville_summary.csv", row.names = FALSE)
# write.csv(trajectory, "eville_trajectory.csv", row.names = FALSE)
# grDevices::pdf("eville_t_test.pdf", width = 7, height = 4.5)
# plot(fit)
# grDevices::dev.off()
