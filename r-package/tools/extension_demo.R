# Run inside the eVille project after devtools::load_all().
# Synthetic illustration; replace x with observations in their actual order.
x <- c(-0.3, -1.1, -0.7, -0.5, -1.4, -0.8, -0.2, -1.0)

prior <- ev_prior("gaussian", location = -4, scale = 1, alternative = "less")
fit <- ev_ttest(x, prior = prior, alpha = 0.05)
print(fit)
print(as.data.frame(fit))

flat <- ev_ttest(x, prior = ev_prior("lebesgue", alternative = "less"))
print(flat)

# Inspect the running p-values. This plot does not change the test.
graphics::matplot(fit$path$n, fit$path[, c("p_ville", "p_eville")],
  type = "l", lty = c(2, 1), lwd = 2,
  col = c("#16324F", "#D97706"), ylim = c(0, 1),
  xlab = "Observations", ylab = "Running anytime-valid p-value")
graphics::lines(flat$path$n, flat$path$p_eville, lty = 3, lwd = 2)
graphics::abline(h = 0.05, lty = 3, col = "grey50")
graphics::legend("topright", c("Ville", "eVille", "Lebesgue eVille"),
  col = c("#16324F", "#D97706", "black"), lty = c(2, 1, 3), bty = "n")
