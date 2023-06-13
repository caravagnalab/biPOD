
#' Plot the fit over the input data.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit' and must have been fitted with model selection.
#' @param plot_type One between "hist", "violin" and "boxplot"
#' @param color Main color of the figure
#'
#' @returns A plot of the Bayes Factor with its significance.
#' @export
plot_mixture_model_omega = function(x, plot_type = "hist", color="maroon") {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")
  if (!("omega_mixture_model" %in% names(x$metadata))) stop("Input must contain a 'omega_mixture_model' field in its metadata!")
  if (!(plot_type %in% c("hist", "violin", "boxplot"))) stop('plot_type must be one of "hist", "violin", "boxplot"')

  omega_samples <- x$metadata$omega_mixture_model %>% dplyr::as_tibble() %>% tidyr::pivot_longer(cols = c("omega"))

  if (plot_type == "hist") {
    p <- ggplot2::ggplot(data=omega_samples, mapping = ggplot2::aes(x=.data$value, y=ggplot2::after_stat(..density..))) +
      ggplot2::geom_histogram(fill = color, col = color, alpha = .6, bins = 100) +
      ggplot2::geom_density(col = color, linewidth = 1.2) +
      biPOD:::my_ggplot_theme() +
      ggplot2::scale_x_continuous(
        name = expression(omega),
        breaks = c(0, 1),
        labels = c("Logistic", "Exponential"),
        limits = c(0, 1)
      ) +
      ggplot2::geom_vline(xintercept = 0.5, color = "darkgray", linetype = "dashed")
  }

  if (plot_type == "violin") {
    p <- ggplot2::ggplot(data=omega_samples, mapping = ggplot2::aes(x=.data$name, y=.data$value, fill=.data$name)) +
      ggplot2::geom_violin() +
      ggplot2::lims(y = c(0,1)) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c(ggplot2::alpha(color, .8))) +
      biPOD:::my_ggplot_theme() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_y_continuous(
        name = expression(omega),
        breaks = seq(0, 1, by = .25),
        labels = c("Logistic", "", "", "", "Exponential"),
        limits = c(0, 1)
      ) +
      ggplot2::scale_x_discrete(
        name = "",
        breaks = NULL
      ) +
      ggplot2::geom_hline(yintercept = 0.5, color = "darkgray", linetype = "dashed")
  }

  if (plot_type == "boxplot") {
    p <- ggplot2::ggplot(data=omega_samples, mapping = ggplot2::aes(x=.data$name, y=.data$value, fill=.data$name)) +
      ggplot2::geom_boxplot() +
      ggplot2::lims(y = c(0,1)) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c(ggplot2::alpha(color, .8))) +
      biPOD:::my_ggplot_theme() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_y_continuous(
        name = expression(omega),
        breaks = seq(0, 1, by = .25),
        labels = c("Logistic", "", "", "", "Exponential"),
        limits = c(0, 1)
      ) +
      ggplot2::scale_x_discrete(
        name = "",
        breaks = NULL
      ) +
      ggplot2::geom_hline(yintercept = 0.5, color = "darkgray", linetype = "dashed")
  }

  p <- p + ggplot2::ggtitle(label = paste0("Odds ratio of ", round(x$metadata$odd, 2), " in favour of ", x$metadata$best_growth))

  return(p)
}
