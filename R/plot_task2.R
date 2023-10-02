#' Plot the posterior distributions of the inferred breakpoints.
#'
#' @param x a bipod object with a 'fit' field
#' @param with_histogram .
#' @param alpha .
#' @param colors colors to use for the different inferred breakpoints
#'
#' @return A ggplot object containing the posterior distributions of the inferred breakpoints.
#' @export
#'
plot_breakpoints_posterior <- function(x, with_histogram = F, alpha = .6, colors = NULL) {
  # plot_breakpoints_posterior function
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("breakpoints_fit" %in% names(x))) stop("Input must contain a 'breakpoints_fit' field")

  n_changepoints <- length(x$metadata$breakpoints)
  par_list <- lapply(1:n_changepoints, function(i) {
    paste0("b[", i, "]")
  }) %>% unlist()

  samples <- get_parameters(x$breakpoints_fit, par_list = par_list)
  samples$value <- samples$value + min(x$counts$time)

  colors = rep('darkgray', length(x$metadata$breakpoints))

  ggplot2::ggplot() +
    ggplot2::geom_density(data = samples, mapping = ggplot2::aes(x = .data$value, fill = .data$parameter, col = .data$parameter), alpha = alpha) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_color_manual(values = colors) +
    my_ggplot_theme()
}
