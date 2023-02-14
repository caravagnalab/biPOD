
#' Plot the input data of the population over time
#'
#' @param x A biPOD object of class `bipod`.
#' @param add_title Boolean, indicating whether the plot should have a title
#' @param log_scale Boolean, indicating whether the plot should have a title
#'
#' @returns A plot. Represents the evolution of the population over time.
#' @export
plot_input <- function(x, add_title = F, log_scale = F) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")

  # If log, transform y scale
  trans <- if(log_scale) "log" else "identity"

  data <- data.frame(x = x$counts$time, y = x$counts$count)

  pop.plot <- ggplot2::ggplot(data, ggplot2::aes(x=.data$x, y=.data$y)) +
    ggplot2::geom_line(col='forestgreen') +
    ggplot2::geom_point(col='forestgreen') +
    ggplot2::scale_y_continuous(trans = trans) +
    my_ggplot_theme()

  if (add_title) {
    pop.plot <- pop.plot + ggplot2::labs(
      title = x$sample,
      x = "time",
      y = "count"
    )
  } else {
    pop.plot <- pop.plot + ggplot2::labs(
      x = "time",
      y = "count"
    )
  }

  pop.plot
}
