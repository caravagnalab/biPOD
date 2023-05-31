
#' Plot the fit over the input data.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param legend_labels Vector containing a name for each unique group in x$counts$group
#' @param legend_title Title for the legend. Default is "group"
#' @param zoom_limits Limits of the x-axis for the zoom of the plot
#' @param sec_axis_breaks Vector containing values in which secondary x axis' breaks will be
#' @param sec_axis_labels Vector containing labels for the secondary x axis
#' @param CI confidence interval for the growth rate to plot
#'
#' @returns A plot of the fit over the input data.
#' @export
plot_fit = function(x,
                    CI = .95,
                    legend_labels = NULL,
                    legend_title = "group",
                    zoom_limits = NULL,
                    sec_axis_breaks = NULL,
                    sec_axis_labels = NULL) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  alpha = 1 - CI
  fitted_data <- get_data_for_plot(x, alpha = alpha)

  # Plot data
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(x$counts, mapping = ggplot2::aes(x=.data$time, y=.data$count)) + #original points
    ggplot2::geom_line(fitted_data, mapping = ggplot2::aes(x=.data$x, y=.data$y), col="black") +
    ggplot2::geom_ribbon(fitted_data, mapping = ggplot2::aes(x=.data$x, y=.data$y, ymin=.data$ylow, ymax=.data$yhigh), fill="black", alpha=.3) +
    biPOD:::my_ggplot_theme()

  # add highlights
  p <- biPOD:::add_shadow_to_plot(x, base_plot = p)

  # add t0 posterior
  p <- biPOD:::add_t0_posterior(base_plot = p, x = x)

  # change legend
  if (!(is.null(legend_labels))) {
    if (!(length(legend_labels) == length(unique(x$counts$group)))) stop("The number of labels for the legend must be equal to the number of groups")
    p <- p + ggplot2::scale_fill_manual(values = get_group_colors(), labels=legend_labels)
  }
  p <- p + ggplot2::guides(fill=ggplot2::guide_legend(title=legend_title, override.aes = list(alpha = 1)))

  # Add zoom
  zoom_limits <- if (is.null(zoom_limits)) c(min(x$counts$time), max(x$counts$time)) else zoom_limits
  p <- p + ggforce::facet_zoom(xlim = zoom_limits)

  # Add secondary axis
  if (xor(is.null(sec_axis_breaks), is.null(sec_axis_labels))) stop("To plot secondary x axis both sec_axis_breaks and sex_axis_labels are needed")
  if (!(is.null(sec_axis_breaks))) {
    p <- p +
      ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~ ., breaks = sec_axis_breaks, labels = sec_axis_labels)) +
      ggplot2::theme(axis.text.x.top = ggplot2::element_text(angle = 90, vjust = 0.5))
  }

  return(p)
}
