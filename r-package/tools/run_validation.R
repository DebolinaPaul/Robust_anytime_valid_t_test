# Run from the package root: source("tools/run_validation.R")
# Runs the full package suite, including shifted and one-sided tests.
if (!file.exists("DESCRIPTION") || !dir.exists("tests/testthat")) {
  stop("Run this script from the R package directory containing DESCRIPTION.")
}
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("First install.packages('testthat').")
}
testthat::test_local(path = ".", reporter = "summary", stop_on_failure = TRUE)
