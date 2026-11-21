# Compare the five simulation priors on one short synthetic sequence.
# Priors and directions are fixed independently of the target observations.
# MC calibration is a finite-sample approximation, with reported diagnostics.
library(eVille)

x <- c(-0.3, 0.8, 1.2, -0.4, 1.1, 0.6)
families <- c("gaussian", "cauchy", "student_t5", "laplace", "horseshoe")

# These are exactly the notebook's default scales. All full-line locations
# are set to 3 for this illustration. No prior is fitted to x.
priors <- stats::setNames(lapply(families, function(family) {
  ev_prior(family, location = 3, alternative = "two.sided")
}), families)
null_t <- ev_null_sample()
calibrators <- lapply(priors, function(prior) {
  ev_calibration(prior, method = "monte_carlo", null_t = null_t)
})
fits <- lapply(families, function(family) {
  ev_ttest(x, prior = priors[[family]], calibration = calibrators[[family]])
})
names(fits) <- families

comparison <- do.call(rbind, lapply(families, function(family) {
  result <- summary(fits[[family]])$results
  data.frame(prior = family, result, row.names = NULL)
}))
print(comparison, row.names = FALSE)

# Plot one prior; all fits share the observations and null sample.
plot(fits$student_t5, main = "Student-t5 prior: Monte Carlo calibration")

# A centered Rouder prior uses the same Cauchy density and default scale.
rouder <- ev_ttest(x, prior = ev_prior("rouder"), null_t = null_t)
flat <- ev_ttest(x, prior = ev_prior("lebesgue"))
print(summary(rouder))
print(summary(flat))

# Density evaluation; the horseshoe density is infinite at its location.
print(ev_prior_density(c(2, 3, 4), priors$horseshoe))
