
#' Plot the evolution of the population over time
#'
#' @param x A biPOD object of class `bipod`.
#' @returns A plot. Represents the evolution of the population over time.
#' @import ggplot2
#' @export
evolution_plot <- function(x) {
  # Check input
  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")

  data <- data.frame(x = x$counts$time, y = x$counts$count)

  pop.plot <- ggplot(data, aes(x=.data$x, y=.data$y)) +
    geom_line(col='forestgreen') +
    geom_point(col='forestgreen') +
    ggtitle(x$sample) +
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
