plot_posterior <- function(x, x_fit, par_name, color = "black") {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  d <- get_parameter(x_fit, par_name)
  p <- ggplot2::ggplot() +
    ggplot2::geom_density(data = d, mapping = ggplot2::aes(x = .data$value), col = "black", fill = color, linewidth = .8, alpha = .6) +
    ggplot2::facet_wrap(~ .data$parameter, labeller = ggplot2::label_parsed) +
    ggplot2::labs(y = "density", x = "") +
    my_ggplot_theme() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = ggplot2::alpha(color, .6)))
  p
}

plot_posteriors <- function(x, x_fit, par_list, with_histogram = F, alpha = .6, colors = NULL) {
  samples <- get_parameters(x_fit, par_list = par_list)

  if (with_histogram) alpha <- 0

  if (is.null(colors)) {
    colors = get_group_colors()
  } else {
    if (length(colors) != length(par_list)) stop("Length of colors should be equal to the length of par_list")
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_density(data = samples, mapping = ggplot2::aes(x = .data$value, fill = .data$parameter, col = .data$parameter), alpha = alpha) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_color_manual(values = colors) +
    my_ggplot_theme()

  if (with_histogram) {
    for (par in par_list) {
      d <- samples %>% dplyr::filter(.data$parameter == par)
      p <- p +
        ggplot2::geom_histogram(data = d, mapping = ggplot2::aes(x = .data$value, y = ggplot2::after_stat(..density..)), alpha = .5, bins = 50)
    }
  }
  p
}
