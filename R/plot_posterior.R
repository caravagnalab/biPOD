
plot_posterior = function(x, x_fit, par_name, color = "black") {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  d <- biPOD:::extract_parameter(x_fit, par_name)
  p <- ggplot2::ggplot() +
    ggplot2::geom_density(data = d, mapping=ggplot2::aes(x=.data$value), col = "black", fill = color, linewidth = .8, alpha = .6) +
    ggplot2::facet_wrap( ~ .data$parameter, labeller = ggplot2::label_parsed) +
    ggplot2::labs( y = 'density', x = '') +
    biPOD:::my_ggplot_theme() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = ggplot2::alpha(color, .6)))
  p
}
