extension_data <- function(n, t) {
  z <- seq_len(n) - (n + 1) / 2
  z / stats::sd(z) + t / sqrt(n)
}

extension_prior_row <- function(row) {
  if (row$prior == 'lebesgue') {
    ev_prior('lebesgue', alternative = row$alternative)
  } else {
    ev_prior(row$prior, location = row$location, scale = row$scale,
             alternative = row$alternative)
  }
}

