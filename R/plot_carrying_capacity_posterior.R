#' Plot posteriors of carrying capacity K
#'
#' @param x a bipod object with a 'fit' field
#' @param add_prior Boolean, indicate whether to plot also the prior distribution
#' @param with_histogram .
#' @param alpha .
#'
#' @return A ggplot object containing the posterior density plots of the growth rates and the prior density plot
#' @export
#'
plot_carrying_capacity_posterior <- function(x, add_prior = F, with_histogram = F, alpha = .6) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")
  if (!(x$metadata$growth_type == "logistic")) stop("Carrying capacity posterior is available only for 'logistic' model!")

  # plot posterior density
  # p <- plot_posterior(x, x$fit, par_name = "K", color = "seagreen")
  p <- plot_posteriors(x, x$fit,
    par_list = c("K"),
    with_histogram = with_histogram, alpha = alpha
  )

  if (add_prior) {
    prior_K <- x$metadata$prior_K
    xmin <- prior_K * .95
    xmax <- prior_K * 10
    xs <- seq(xmin, xmax, length = 500)
    ys <- stats::dnorm(xs, mean = prior_K, sd = prior_K)
    prior_data <- dplyr::tibble(x = xs, y = ys)
    p <- p + ggplot2::geom_line(
      data = prior_data,
      ggplot2::aes(x = .data$x, y = .data$y),
      col = "indianred3",
      size = .8
    )
  }

  # Add style
  p <- p +
    ggplot2::labs(y = "density", x = "value") +
    # ggplot2::coord_cartesian(xlim=xlims) +
    my_ggplot_theme()

  p
}
