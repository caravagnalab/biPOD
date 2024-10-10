#' Plot Posterior Distributions of Inferred Breakpoints
#'
#' Generates a plot showing the posterior distributions of the inferred breakpoints from a fitted model.
#' The plot can include density estimates and optional histograms for better visualization.
#'
#' @param x A `bipod` object that contains a 'breakpoints_fit' field.
#' @param with_histogram A logical value indicating whether to overlay histograms of the breakpoint samples on the density plots. (default is FALSE)
#' @param alpha A numeric value between 0 and 1 specifying the transparency level of the density plot. (default is 0.6)
#' @param colors A character vector specifying colors to use for the different inferred breakpoints. If `NULL`, the default color 'darkgray' is used for all breakpoints.
#'
#' @return A `ggplot2` object displaying the posterior distributions of the inferred breakpoints. The plot shows density estimates and optionally histograms for the breakpoints.
#' @export
plot_breakpoints_posterior <- function(x, with_histogram = F, alpha = .6, colors = NULL) {
  # plot_breakpoints_posterior function
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("breakpoints_fit" %in% names(x))) stop("Input must contain a 'breakpoints_fit' field")

  n_changepoints <- length(x$metadata$breakpoints)
  par_list <- lapply(1:n_changepoints, function(i) {
    paste0("b[", i, "]")
  }) %>% unlist()

  samples <- get_parameters(x$breakpoints_fit, par_list = par_list, variational = F)
  if (!(is.null(x$breakpoints_fit$normalization_pars))) {
    samples$value <- samples$value * x$breakpoints_fit$normalization_pars$sd + x$breakpoints_fit$normalization_pars$mean
  }

  if (is.null(colors)) { colors = rep('darkgray', length(x$metadata$breakpoints)) }

  ggplot2::ggplot() +
    ggplot2::geom_density(data = samples, mapping = ggplot2::aes(x = .data$value, fill = .data$parameter, col = .data$parameter), alpha = alpha) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_color_manual(values = colors) +
    my_ggplot_theme()
}
