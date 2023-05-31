
#' Plot posteriors of t0
#'
#' @param x a bipod object with a 'fit' field
#' @param add_prior Boolean, indicate whether to plot also the prior distribution
#'
#' @return A ggplot object containing the posterior density plots of the growth rates and the prior density plot
#' @export
#'
plot_t0_posterior = function(x, add_prior = F) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  t0_lower_bound <- x$metadata$t0_lower_bound

  # Plot model with t0
  if (t0_lower_bound == x$counts$time[1]) {
    stop("t0 has been set manually as the first time step. Posterior distribution not available.")
  } else {
    p <- biPOD:::plot_posterior(x, x$fit, "t0", "darkorange")
  }

  if (add_prior) {
    xmin <- x$metadata$t0_lower_bound - 0.1
    xmax <- min(x$counts$time) + 0.1
    xs <- seq(xmin, xmax, length=500)
    ys <- stats::dunif(xs, x$metadata$t0_lower_bound, min(x$counts$time))
    prior_data = dplyr::tibble(x=xs, y=ys)
    p <- p + ggplot2::geom_line(
      data = prior_data,
      ggplot2::aes(x=.data$x, y=.data$y),
      col = "indianred3",
      size = .8)
  }

  # Add style
  p <- p +
    ggplot2::labs( y = 'density', x = '') +
    biPOD:::my_ggplot_theme() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = ggplot2::alpha("darkorange", .8)))

  p
}
