#' Plot the input data of the population over time
#'
#' @param x A bipod object of class `bipod`.
#' @param log_scale Boolean, indicating whether the plot should have a title
#' @param add_highlights Boolean, indicating whether the groups should be highlighted
#'
#' @returns A plot. Represents the evolution of the population over time.
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
    biPOD:::my_ggplot_theme()

  if (add_highlights) {
    p <- biPOD:::add_shadow_to_plot(x, base_plot = p)
  }

  # add style
  p <- p + ggplot2::labs(
    x = "time",
    y = y_label
  )

  p
}
