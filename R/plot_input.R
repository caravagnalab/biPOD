
#' Plot Population Counts Over Time
#'
#' Generates a plot of population counts over time from a `bipod` object. This function allows customization of the plot, including the option to display counts on a logarithmic scale and to highlight different groups or categories within the data.
#' The plot visualizes the evolution of population counts at each time point, helping to understand trends and fluctuations over time.
#'
#' @param x A `bipod` object containing population count data. The object must include a `counts` field with columns `time` (time points) and `count` (population counts).
#' @param log_scale A logical value indicating whether to apply a logarithmic scale to the y-axis. If `TRUE`, the y-axis will be transformed to a logarithmic scale, which can help visualize data with large differences in magnitude. If `FALSE` (default), the y-axis will use a linear scale.
#' @param add_highlights A logical value indicating whether to highlight different groups within the population data. If `TRUE`, the plot will include additional visual elements to differentiate between groups or categories. If `FALSE` (default), no highlighting will be applied.
#'
#' @return A `ggplot2` object representing the evolution of the population counts over time. The plot will include lines and points for the population counts and can be customized based on the provided arguments.
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
