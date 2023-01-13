
#' Plot the evolution of the population over time
#'
#' @param x A biPOD object of class `bipod`.
#'
#' @returns A plot. Represents the evolution of the population over time.
#' @export
evolution_plot <- function(x) {
  # Check input
  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")

  data <- data.frame(x = x$counts$time, y = x$counts$count)

  pop.plot <- ggplot2::ggplot(data, ggplot2::aes(x=.data$x, y=.data$y)) +
    ggplot2::geom_line(col='forestgreen') +
    ggplot2::geom_point(col='forestgreen') +
    ggplot2::labs(
      title = x$sample
    )
    my_ggplot_theme()

  pop.plot
}

my_ggplot_theme = function(cex_opt = 1)
{
  ggplot2::theme_light(base_size = 10 * cex_opt) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(.3 * cex_opt, "cm"),
      panel.background = ggplot2::element_rect(fill = 'white')
    )
}
