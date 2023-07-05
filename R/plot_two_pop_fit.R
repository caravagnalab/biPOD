#' Plot the fit over the input data.
#'
#' @param x A bipod object of class `bipod`. Must contains 'fit'
#' @param CI confidence interval for the growth rate to plot
#'
#' @returns A plot of the fit over the input data.
#' @export
plot_two_pop_fit <- function(x, CI = .95) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("two_pop_fit" %in% names(x))) stop("Input must contain a 'two_pop_fit' field")

  alpha <- 1 - CI
  fitted_data <- get_data_for_two_pop_plot(x, alpha = alpha)

  # Plot data
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(x$counts, mapping = ggplot2::aes(x = .data$time, y = .data$count)) + # original points
    ggplot2::geom_line(fitted_data %>% dplyr::filter(.data$group == "resistant"), mapping = ggplot2::aes(x = .data$x, y = .data$y), col = "forestgreen") +
    ggplot2::geom_ribbon(fitted_data %>% dplyr::filter(.data$group == "resistant"), mapping = ggplot2::aes(x = .data$x, y = .data$y, ymin = .data$ylow, ymax = .data$yhigh), fill = "forestgreen", alpha = .3) +
    ggplot2::geom_line(fitted_data %>% dplyr::filter(.data$group == "sensible"), mapping = ggplot2::aes(x = .data$x, y = .data$y), col = "indianred") +
    ggplot2::geom_ribbon(fitted_data %>% dplyr::filter(.data$group == "sensible"), mapping = ggplot2::aes(x = .data$x, y = .data$y, ymin = .data$ylow, ymax = .data$yhigh), fill = "indianred", alpha = .3) +
    my_ggplot_theme()

  p <- add_posterior(base_plot = p, param = "t0_r", x = x, x_fit = x$two_pop_fit, color = "forestgreen")
  p <- add_posterior(base_plot = p, param = "t_end", x = x, x_fit = x$two_pop_fit, color = "indianred")

  return(p)
}
