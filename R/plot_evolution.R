
#' Plot the evolution of the population over time
#'
#' @param x A biPOD object of class `bipod`.
#' @param add_title Boolean, indicating whether the plot should have a title
#'
#' @returns A plot. Represents the evolution of the population over time.
#' @export
plot_evolution <- function(x, add_title = F) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")

  data <- data.frame(x = x$counts$time, y = x$counts$count)

  pop.plot <- ggplot2::ggplot(data, ggplot2::aes(x=.data$x, y=.data$y)) +
    ggplot2::geom_line(col='forestgreen') +
    ggplot2::geom_point(col='forestgreen') +
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
