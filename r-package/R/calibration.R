#' Two-sided Lebesgue calibration starting at n = 2
#'
#' @param c Numeric vector of nonnegative evidence thresholds.
#'   Infinite thresholds are allowed.
#'
#' @return A numeric vector with values between 0 and 1 (inclusive).
#' @export
#'
#' @examples
#' g_lebesgue(c(0, sqrt(pi), 10, Inf))
g_lebesgue <- function(c) {
  if (!is.numeric(c) || anyNA(c) || any(c < 0)) {
    stop(
      "`c` must contain nonnegative numbers without missing values.",
      call. = FALSE
    )
  }

  result <- rep(1, length(c))
  above <- c > sqrt(pi)

  q <- sqrt(pi) / c[above]

  result[above] <- (2 / pi) * (
    asin(sqrt(q)) + sqrt(q * (1 - q))
  )

  result
}
