#' Plot Posterior Density for a Specific Parameter
#'
#' Creates a density plot of the posterior distribution for a specified parameter from a fitted model within a `bipod` object.
#' The plot helps visualize the distribution of the parameter values after fitting the model, indicating the uncertainty in the parameter estimate.
#' This function is particularly useful for exploring the posterior distribution of individual model parameters in Bayesian analysis.
#'
#' @param x A `bipod` object of class `bipod`. It must contain a 'fit' field, which includes the fitted model results.
#'          This object should have been fitted using model selection techniques.
#' @param x_fit The fit object within `x` that contains the desired parameter. It should include the fitted parameter values and associated diagnostics.
#' @param par_name A character string specifying the name of the parameter whose posterior distribution is to be plotted.
#'                 This should match one of the parameter names available in the model fit.
#' @param color A character string specifying the color to use for the density plot. The default is `'black'`.
#'              This color will be used for filling the density plot, with a specified level of transparency.
#'
#' @return A `ggplot2` object representing the density plot of the specified parameter. The plot displays the posterior distribution
#'         and helps assess the uncertainty of the parameter estimate.
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
#' Generates density plots for the posterior distributions of multiple parameters from a fitted model within a `bipod` object.
#' This function can optionally overlay histograms of the posterior samples on top of the density plots for a better visual comparison of the distribution.
#' This is useful when you want to examine the posterior distributions of several parameters simultaneously.
#'
#' @param x A `bipod` object of class `bipod`. It must contain a 'fit' field with the fitted model and associated results.
#'          This object should have undergone model selection and fitting.
#' @param x_fit The fit object within `x` that contains the desired parameters. It holds the fitted values for the parameters to be plotted.
#' @param par_list A character vector containing the names of the parameters whose posterior distributions are to be plotted.
#'                 These should match the parameter names available in the model fit.
#' @param with_histogram A logical value indicating whether to overlay histograms on the density plots.
#'                       If `TRUE`, histograms will be displayed alongside the density plots to give a clearer visualization of the sample distribution.
#'                       The default is `FALSE`.
#' @param alpha A numeric value between 0 and 1 specifying the transparency level of the density plots.
#'              A lower value (e.g., 0.2) will make the density plot more transparent, while a higher value (e.g., 0.8) will make it more opaque. The default is `0.6`.
#' @param colors A character vector specifying the colors to use for the density plots of each parameter.
#'               If `NULL` (default), a predefined color scheme will be used.
#'               The length of this vector should match the number of parameters in `par_list`.
#'
#' @return A `ggplot2` object displaying the density plots of the specified parameters.
#'         If `with_histogram` is `TRUE`, histograms will be overlaid on the density plots.
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
