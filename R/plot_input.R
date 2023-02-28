
#' Plot the input data of the population over time
#'
#' @param x A biPOD object of class `bipod`.
#' @param add_title Boolean, indicating whether the plot should have a title
#' @param log_scale Boolean, indicating whether the plot should have a title
#' @param add_highlights Boolean, indicating whether the groups should be highlighted
#'
#' @returns A plot. Represents the evolution of the population over time.
#' @export
plot_input <- function(x, add_title = F, log_scale = F, add_highlights = F) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")

  # If log, transform y scale
  trans <- if(log_scale) "log" else "identity"
  y_label <- if(log_scale) "log count" else "count"

  data <- data.frame(x = x$counts$time, y = x$counts$count)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data, mapping=ggplot2::aes(x=.data$x, y=.data$y), col='black') +
    ggplot2::geom_point(data, mapping=ggplot2::aes(x=.data$x, y=.data$y), col='black') +
    ggplot2::scale_y_continuous(trans = trans) +
    my_ggplot_theme()

  if (add_highlights) {
    p <- add_shadow_to_plot(x, base_plot = p)
  }

  if (add_title) {
    p <- p + ggplot2::labs(
      title = paste("Sample:", x$sample),
      x = "time",
      y = y_label
    )
  } else {
    p <- p + ggplot2::labs(
      x = "time",
      y = y_label
    )
  }

  p
}
