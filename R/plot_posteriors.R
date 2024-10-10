#' Plot Posterior Density for a Specific Parameter
#'
#' Creates a density plot of the posterior distribution for a specified parameter from a fitted model within a `bipod` object.
#'
#' @param x A `bipod` object of class `bipod`. Must contain a 'fit' field and should have been fitted with model selection.
#' @param x_fit The fit object within `x` that contains the desired parameter.
#' @param par_name The name of the parameter whose posterior distribution is to be plotted.
#' @param color A character string specifying the color to use for the density plot. (default is 'black')
#'
#' @return A `ggplot2` object showing the density plot of the specified parameter.
#' @export
plot_posterior <- function(x, x_fit, par_name, color = "black") {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!(par_name %in% x_fit$parameters)) stop(paste0("Input par_name is not an available parameters for the given input fit.\nAvailable parameters are : ", paste(x_fit$parameters, collapse = ", ")))
  d <- get_parameter(x_fit, par_name, x$metadata$sampling == "variational")
  p <- ggplot2::ggplot() +
    ggplot2::geom_density(data = d, mapping = ggplot2::aes(x = .data$value), col = "black", fill = color, linewidth = .8, alpha = .6) +
    ggplot2::facet_wrap(~ .data$parameter, labeller = ggplot2::label_parsed) +
    ggplot2::labs(y = "density", x = "") +
    my_ggplot_theme() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = ggplot2::alpha(color, .6)))
  p
}

#' Plot Posterior Densities for Multiple Parameters
#'
#' Generates density plots for the posterior distributions of a set of parameters from a fitted model within a `bipod` object.
#'  Optionally overlays histograms on the density plots.
#'
#' @param x A `bipod` object of class `bipod`. Must contain a 'fit' field and should have been fitted with model selection.
#' @param x_fit The fit object within `x` that contains the desired parameters.
#' @param par_list A character vector containing the names of the parameters whose posterior distributions are to be plotted.
#' @param with_histogram A logical value indicating whether to overlay histograms of the posterior samples on the density plots. (default is FALSE)
#' @param alpha A numeric value between 0 and 1 specifying the transparency level of the density plots. (default is .6)
#' @param colors A character vector specifying the colors to use for the density plots of each parameter.
#'  If `NULL`, a default color scheme will be used. (default is NULL)
#'
#' @return A `ggplot2` object showing the density plots of the specified parameters. If `with_histogram` is `TRUE`, histograms will be overlaid on the density plots.
#' @export
plot_posteriors <- function(x, x_fit, par_list, with_histogram = F, alpha = .6, colors = NULL) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")

  samples <- get_parameters(x_fit, par_list = par_list, variational = x$metadata$sampling == "variational")

  if (with_histogram) alpha <- 0

  if (is.null(colors)) {
    colors <- get_group_colors()
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
