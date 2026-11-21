#' Summarise an anytime-valid t-test
#'
#' @param object An object returned by ev_ttest().
#' @param ... Reserved for future extensions.
#' @return A summary.ev_ttest object containing the final running p-values,
#'   first rejection times, prior settings, and calibration diagnostics.
#' @details The summary uses the stored results without recomputing evidence
#'   or calibration. An NA rejection time means that no rejection was observed,
#'   except for Ville with Lebesgue measure, where the method is unavailable.
#'   Prior location and scale refer to the family before one-sided conditioning.
#'   The calibration integration gap does not include floating-point error in
#'   the evidence calculation. In Monte Carlo mode, the summary also reports
#'   the null sample size and empirical calibration standard error. These
#'   diagnostics do not establish exact anytime validity for finite Monte
#'   Carlo calibration.
#' @export
summary.ev_ttest <- function(object, ...) {
  path <- object$path
  last <- nrow(path)
  available <- c(TRUE, isTRUE(object$prior$proper))
  first <- unname(object$stopping_times[c("eVille", "Ville")])
  p <- c(path$p_eville[last], path$p_ville[last])
  p[!available] <- NA_real_
  first[!available] <- NA_integer_
  decision <- ifelse(!available, "Unavailable",
                     ifelse(is.na(first), "No rejection observed", "Rejected"))
  results <- data.frame(
    method = c("eVille", "Ville"),
    final_running_p = p,
    first_rejection = first,
    decision = decision
  )
  structure(list(
    call = object$call,
    prior = object$prior,
    mu = object$mu,
    alpha = object$alpha,
    n = path$n[last],
    n0 = path$n[1L],
    results = results,
    max_calibration_gap = max(path$calibration_gap),
    calibration_method = object$calibration$method,
    mc_samples = if (identical(object$calibration$method, "monte_carlo")) {
      object$calibration$B
    } else 0L,
    final_mcse = if ("mcse" %in% names(path)) path$mcse[last] else NA_real_,
    max_mc_gap = if ("mc_gap" %in% names(path)) max(path$mc_gap) else NA_real_,
    calibration_ambiguous_n = object$calibration_ambiguous_n
  ), class = "summary.ev_ttest")
}

#' Print an anytime-valid t-test summary
#'
#' @param x A summary returned by summary() for an ev_ttest object.
#' @param digits Number of significant digits displayed.
#' @param ... Reserved for future extensions.
#' @return The input, invisibly.
#' @export
print.summary.ev_ttest <- function(x, digits = 4, ...) {
  mc <- identical(x$calibration_method, "monte_carlo")
  cat(if (mc) "One-sample t-test (Monte Carlo eVille calibration)\n" else
    "Anytime-valid one-sample t-test\n")
  cat("Observations:", x$n, "| Monitoring starts at:", x$n0, "\n")
  cat("Calibration:", x$calibration_method, "\n")
  cat("Alternative:", x$prior$alternative,
      "| Null boundary:", x$mu, "| Alpha:", x$alpha, "\n")
  cat("Prior:", x$prior$family)
  if (x$prior$proper) {
    cat(" | location:", x$prior$location, "| scale:", x$prior$scale)
    if (x$prior$alternative != "two.sided") {
      cat(" (before one-sided conditioning)")
    }
  }
  cat("\n\n")
  display <- x$results
  display$final_running_p <- format(display$final_running_p,
                                    digits = digits, trim = TRUE)
  display$first_rejection <- ifelse(
    is.na(display$first_rejection), "Not reached",
    as.character(display$first_rejection)
  )
  unavailable <- display$decision == "Unavailable"
  display$final_running_p[unavailable] <- "--"
  display$first_rejection[unavailable] <- "--"
  print(display, row.names = FALSE, right = FALSE)
  cat("\nMaximum numerical calibration gap:",
      format(x$max_calibration_gap, digits = digits), "\n")
  if (mc) {
    cat("Null Monte Carlo draws:", x$mc_samples, "\n")
    cat("Final empirical calibration MC standard error (before Ville cap):",
        format(x$final_mcse, digits = digits), "\n")
    cat("Maximum raw empirical numerical gap:",
        format(x$max_mc_gap, digits = digits), "\n")
    cat("Finite Monte Carlo approximation; MC error is separate from the numerical gap.\n")
  }
  if (length(x$calibration_ambiguous_n)) {
    cat("Calibration bracket straddles alpha at n =",
        paste(x$calibration_ambiguous_n, collapse = ", "), "\n")
    cat("Decisions use the upper calibrated p-value; refine the calibration",
        "to resolve these times.\n")
  }
  invisible(x)
}

#' Plot running anytime-valid p-values
#'
#' @param x An object returned by ev_ttest().
#' @param y Unused; included for compatibility with plot().
#' @param ... Additional graphical arguments, such as xlim or cex.axis,
#'   passed to graphics::plot(). The type and log arguments are not supported:
#'   this method draws step paths on linear axes.
#' @param method Which procedure to display: "both", "eVille", or "Ville".
#' @param xlab,ylab Axis labels.
#' @param main Plot title. The default identifies Monte Carlo calibration when used.
#' @param ylim Vertical axis limits, default c(0, 1).
#' @param legend_position A legend position such as "topright", or NULL
#'   to omit the legend.
#' @param mark_rejections Whether to mark each observed first rejection.
#' @return The input, invisibly. A plot is drawn on the current graphics device.
#' @details The plot uses the stored running p-values. The horizontal line
#'   shows the significance level used in the fitted test. With Lebesgue
#'   measure, "both" displays only eVille; requesting "Ville" gives an error.
#' @importFrom graphics plot
#' @export
plot.ev_ttest <- function(x, y = NULL, ...,
                         method = c("both", "eVille", "Ville"),
                         xlab = "Number of observations",
                         ylab = "Running anytime-valid p-value",
                         main = "Anytime-valid one-sample t-test", ylim = c(0, 1),
                         legend_position = "topright",
                         mark_rejections = TRUE) {
  if (any(c("type", "log") %in% names(list(...)))) {
    stop("This method uses step paths on linear axes; omit type and log.",
         call. = FALSE)
  }
  method <- match.arg(method)
  mc <- identical(x$calibration$method, "monte_carlo")
  if (mc && missing(main)) main <- "One-sample t-test: Monte Carlo eVille calibration"
  if (mc && missing(ylab)) ylab <- "Running p-value (Monte Carlo eVille)"
  selected <- if (method == "both") c("eVille", "Ville") else method
  if (!isTRUE(x$prior$proper)) {
    if (method == "Ville") {
      stop("Ville calibration is unavailable for Lebesgue measure.",
           call. = FALSE)
    }
    selected <- "eVille"
  }
  path <- x$path
  columns <- c(eVille = "p_eville", Ville = "p_ville")
  colours <- c(eVille = "#0072B2", Ville = "#D55E00")
  styles <- c(eVille = 1, Ville = 2)
  graphics::plot(path$n, path[[columns[selected[1L]]]], type = "n",
                 xlab = xlab, ylab = ylab, main = main, ylim = ylim, ...)
  graphics::abline(h = x$alpha, lty = 3, col = "grey40")
  for (name in selected) {
    values <- path[[columns[name]]]
    graphics::lines(path$n, values, type = "s", lwd = 2,
                    col = colours[name], lty = styles[name])
    if (nrow(path) == 1L) {
      graphics::points(path$n, values, pch = 16, col = colours[name])
    }
    first <- unname(x$stopping_times[name])
    if (isTRUE(mark_rejections) && !is.na(first)) {
      j <- match(first, path$n)
      graphics::points(path$n[j], values[j], pch = 16,
                       col = colours[name])
    }
  }
  if (!is.null(legend_position)) {
    graphics::legend(legend_position,
                     legend = c(ifelse(mc & selected == "eVille", "eVille (MC)", selected),
                                paste0("alpha = ", x$alpha)),
                     col = c(unname(colours[selected]), "grey40"),
                     lty = c(unname(styles[selected]), 3),
                     lwd = c(rep(2, length(selected)), 1), bty = "n")
  }
  invisible(x)
}
