
#' Plot the posterior distributions of the inferred breakpoints.
#'
#' @param x a bipod object with a 'fit' field
#'
#' @return A ggplot object containing the posterior distributions of the inferred breakpoints.
#' @export
#'
plot_breakpoints_posterior = function(x) {
  # plot_breakpoints_posterior function
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("breakpoints_fit" %in% names(x))) stop("Input must contain a 'breakpoints_fit' field")

  breakpoints_names <- names(x$breakpoints_fit)[grepl("changing_times\\[", names(x$breakpoints_fit))]
  breakpoints_sample <- biPOD:::extract_parameters(x$breakpoints_fit, par_list = breakpoints_names)

  p <- ggplot2::ggplot() +
    ggplot2::geom_density(data=breakpoints_sample, mapping=ggplot2::aes(x=.data$value, fill=.data$parameter, col=.data$parameter), alpha=.7) +
    ggplot2::scale_fill_manual(values = biPOD:::get_group_colors()[1:length(breakpoints_names)]) +
    ggplot2::scale_color_manual(values = biPOD:::get_group_colors()[1:length(breakpoints_names)]) +
    biPOD:::my_ggplot_theme()
  p
}
