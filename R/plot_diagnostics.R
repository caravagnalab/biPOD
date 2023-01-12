
#' Plot the posterior predictive checks for counts data.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param ptype A string representing the type of visualization.
#'
#' `intervals` vertical bars with points indicating generated values medians and darker points indicating observed values
#' `ribbon` ribbon of connected intervals with a line through the median of generated values and a darker line connecting observed values
#'
#' @param prob,prob_outer Values between 0 and 1 indicating the desired probability mass to include in the inner and outer intervals. The defaults are prob=0.5 and prob_outer=0.9.
#'
#' @returns A plot. Represents the posterior predictive checks.
#' @import ggplot2
#' @importFrom rstan extract
#' @importFrom bayesplot ppc_intervals ppc_ribbon
#' @export
ppc_plot = function(x, ptype, prob = 0.5, prob_outer = 0.9) {
  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  assertthat::assert_that("fit" %in% names(x), msg = "Input must contain a 'fit' field")
  assertthat::assert_that(ptype %in% c("intervals", "ribbon"), msg = "ptype must be one of: intervals, ribbon")
  assertthat::assert_that(prob >= 0 && prob < 1, msg = "prob should be between 0 and 1 (with 1 excluded)")
  assertthat::assert_that(prob_outer >= 0 && prob_outer <= 1, msg = "prob_outer should be between 0 and 1")
  assertthat::assert_that(length(grep("N_rep", names(x$fit))) >= 1, msg = "ppc_plot not support ppc_plot does not support the type of model used")

  y <- x$counts$count[2:length(x$counts$count)]
  yrep <- extract(x$fit, pars="N_rep")$N_rep

  if (ptype == "intervals") {
    p <- ppc_intervals(y, yrep, prob = prob, prob_outer = prob_outer)
  } else {
    p <- ppc_ribbon(y, yrep, prob = prob, prob_outer = prob_outer, y_draw = "points")
  }

  p <- p +
    xlab("Time") +
    ylab("Count") +
    my_ggplot_theme()
  p
}
