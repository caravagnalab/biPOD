
#' Plot the evolution of the population over time
#'
#' @param x A biPOD object of class `bipod`.
#'
#' @returns A plot. Represents the evolution of the population over time.
#' @export
evolution_plot <- function(x) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")

  data <- data.frame(x = x$counts$time, y = x$counts$count)

  pop.plot <- ggplot2::ggplot(data, ggplot2::aes(x=.data$x, y=.data$y)) +
    ggplot2::geom_line(col='forestgreen') +
    ggplot2::geom_point(col='forestgreen') +
    ggplot2::labs(
      title = x$sample,
      x = "Time",
      y = "Count"
    ) +
    my_ggplot_theme()

  pop.plot
}
