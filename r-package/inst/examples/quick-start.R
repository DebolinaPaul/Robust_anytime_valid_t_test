# Anytime-valid one-sample t-test: synthetic worked example.
# Run after installing eVille. Observations remain in their stated order.
library(eVille)

x <- c(-0.3, -1.1, -0.7, -0.5, -1.4, -0.8, -0.2, -1.0)

# Prior specified in advance on theta = (mu - mu0) / sigma.
# This is N(-4, 1), conditioned on theta < 0.
prior <- ev_prior("gaussian", location = -4, scale = 1,
                  alternative = "less")

# H0: mu >= 0 versus H1: mu < 0.
fit <- ev_ttest(x, prior = prior, mu = 0, alpha = 0.05)
print(summary(fit))
plot(fit, main = "Anytime-valid one-sample t-test")

# Numeric summary table and complete sequential trajectory.
summary_table <- summary(fit)$results
trajectory <- as.data.frame(fit)
print(summary_table)

# Separately calibrated half-Lebesgue benchmark, with the same direction.
flat <- ev_ttest(x, prior = ev_prior("lebesgue", alternative = "less"),
                 mu = 0, alpha = 0.05)
print(summary(flat))

# Optional exports: uncomment the following lines if wanted.
# write.csv(summary_table, "eville_summary.csv", row.names = FALSE)
# write.csv(trajectory, "eville_trajectory.csv", row.names = FALSE)
# grDevices::pdf("eville_t_test.pdf", width = 7, height = 4.5)
# plot(fit, main = "Anytime-valid one-sample t-test")
# grDevices::dev.off()
