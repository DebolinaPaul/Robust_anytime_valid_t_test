# Optional diagnostic; deliberately not a pass/fail unit test.
# Load the package with devtools::load_all() before running this function.
null_diagnostic <- function(B = 1000L, n_max = 100L, alpha = 0.05,
                            scale = 30, seed = 20261121L) {
  if (B < 1 || B != as.integer(B) || n_max < 2 ||
      n_max != as.integer(n_max) || alpha <= 0 || alpha >= 1) {
    stop("Require positive integer B, integer n_max >= 2, and 0 < alpha < 1.")
  }
  set.seed(seed)
  rejected <- replicate(B, {
    x <- stats::rnorm(n_max)
    p <- ev_path(x, "gaussian", scale = scale)
    c(Ville = any(p$p_ville <= alpha), eVille = any(p$p_eville <= alpha))
  })
  result <- do.call(rbind, lapply(seq_len(nrow(rejected)), function(j) {
    count <- sum(rejected[j, ])
    interval <- stats::binom.test(count, B)$conf.int
    data.frame(method = rownames(rejected)[j], rejections = count,
               simulations = B, estimate = count / B,
               lower = interval[1], upper = interval[2], alpha = alpha)
  }))
  rownames(result) <- NULL
  result
}
