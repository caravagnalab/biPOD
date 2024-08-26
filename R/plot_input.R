
#' Plot Population Counts Over Time
#'
#' Generates a plot of population counts over time from a `bipod` object.
#' The plot can be customized to display counts on a logarithmic scale and to highlight different groups if specified.
#'
#' @param x A `bipod` object.
#' @param log_scale A logical value indicating whether to apply a logarithmic scale to the y-axis.
#'  If `TRUE`, the y-axis will be transformed to a log scale; otherwise, it will use a linear scale. (default is FALSE)
#' @param add_highlights A logical value indicating whether to highlight different groups in the plot.
#'  If `TRUE`, additional visual elements will be added to the plot to distinguish groups. (defaults if FALSE)
#'
#' @return A `ggplot2` object representing the evolution of the population counts over time.
#' @export
plot_input <- function(x, log_scale = F, add_highlights = F) {
  # Check input
  if (!(inherits(x, "bipod"))) cli::cli_abort(c("{.var x} must be a 'bipod' object", "x" = "You've supplied a {.cls {class(x)}} object"))

  # If log, transform y scale
  trans <- if (log_scale) "log" else "identity"
  y_label <- if (log_scale) "log count" else "count"

  data <- dplyr::tibble(xs = x$counts$time, ys = x$counts$count)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data, mapping = ggplot2::aes(x = .data$xs, y = .data$ys), col = "black") +
    ggplot2::geom_point(data, mapping = ggplot2::aes(x = .data$xs, y = .data$ys), col = "black") +
    ggplot2::scale_y_continuous(trans = trans) +
    my_ggplot_theme()

  if (add_highlights) {
    p <- add_shadow_to_plot(x, base_plot = p, colors = NULL)
  }

  # add style
  p <- p + ggplot2::labs(
    x = "time",
    y = y_label
  )

  p
}
