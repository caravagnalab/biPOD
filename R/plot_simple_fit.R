#' Plot a simple fit over the input data.
#'
#' @param x A bipod object of class `bipod`. Must contains 'fit'
#' @param legend_labels Vector containing a name for each unique group in x$counts$group
#' @param legend_title Title for the legend. Default is "group"
#' @param full_process Boolean, indicating whether to plot the process starting from t0 or not.
#' @param CI confidence interval for the growth rate to plot
#' @param shadows_colors colors to use for the different time windows
#' @param t0_posterior_color .
#'
#' @returns A plot of the fit over the input data.
#' @export
plot_simple_fit <- function(x,
                            CI = .95,
                            full_process = F,
                            legend_labels = NULL,
                            legend_title = "group",
                            t0_posterior_color = "darkorange",
                            shadows_colors = NULL
                            ) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  growth_type <- x$metadata$growth_type
  alpha <- 1 - CI

  fitted_data <- get_data_for_plot(x, alpha = alpha)

  xmin <- min(x$counts$time)
  xmax <- max(x$counts$time)
  if (!(full_process)) {
    fitted_data <- fitted_data %>%
      dplyr::filter(.data$x <= xmax) %>%
      dplyr::filter(.data$x >= xmin)
  }

  # Plot data
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(x$counts, mapping = ggplot2::aes(x = .data$time, y = .data$count)) + # original points
    ggplot2::geom_line(fitted_data, mapping = ggplot2::aes(x = .data$x, y = .data$y), col = "black") +
    ggplot2::geom_ribbon(fitted_data, mapping = ggplot2::aes(x = .data$x, y = .data$y, ymin = .data$ylow, ymax = .data$yhigh), fill = "black", alpha = .3) +
    my_ggplot_theme()

  # add highlights
  p <- add_shadow_to_plot(x, base_plot = p, colors = shadows_colors)

  # add t0 posterior
  if (full_process) p <- add_t0_posterior(base_plot = p, x = x, color=t0_posterior_color)

  # change legend
  if (!(is.null(legend_labels))) {
    if (!(length(legend_labels) == length(unique(x$counts$group)))) stop("The number of labels for the legend must be equal to the number of groups")
    p <- p + ggplot2::scale_fill_manual(values = get_group_colors(), labels = legend_labels)
  }
  p <- p + ggplot2::guides(fill = ggplot2::guide_legend(title = legend_title, override.aes = list(alpha = 1)))

  # Remove everyithing except process eventually
  if (!(full_process)) p <- p + ggplot2::scale_x_continuous(limits = c(xmin, xmax))

  return(p)
}
