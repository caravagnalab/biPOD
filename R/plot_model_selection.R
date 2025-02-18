#' Plot Bayes Factor Over Input Data
#'
#' Generates a plot of the Bayes Factor (BF) based on the fitted model from a `bipod` object. The plot can include categories for interpreting the Bayes Factor according to Jeffreys' scale of evidence.
#'
#' @param x A `bipod` object that contains the results of a fitted model, including Bayes Factor and best growth model metadata.
#' @param with_categories A logical value indicating whether to include Bayes Factor significance categories based on Jeffreys' scale.
#'  If `TRUE`, the plot will use categories to interpret the Bayes Factor. (default is FALSE)
#'
#' @return A `ggplot2` object displaying the Bayes Factor with its significance, optionally categorized.
#'
#' @export
plot_bayes_factor <- function(x, with_categories = F) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")
  if (!("bayes_factor" %in% names(x$metadata))) stop("Input must contain a 'bayes_factor' field in its metadata!")
  if (!("best_growth" %in% names(x$metadata))) stop("Input must contain a 'best_growth' field in its metadata!")

  logBF <- x$metadata$bayes_factor %>%
    as.numeric() %>%
    log10()
  if (with_categories) {
    logBF <- min(3, logBF)
  } else {
    logBF <- min(10, logBF)
  }

  if (x$metadata$best_growth == "Exponential") logBF <- logBF * (-1)

  if (with_categories) {
    p <- empty_bayes_factor_plot_with_categories() +
      ggplot2::geom_point(data = NULL, mapping = ggplot2::aes(x = logBF, y = 0), size = 2)
  } else {
    p <- empty_bayes_factor_plot() +
      ggplot2::geom_point(data = NULL, mapping = ggplot2::aes(x = logBF, y = 0), size = 2)
  }

  p
}

empty_bayes_factor_plot_with_categories <- function() {
  x_from <- c(-3, -2, -3 / 2, -1, -1 / 2, 1 / 2, 1, 3 / 2, 2)
  x_to <- c(-2, -3 / 2, -1, -1 / 2, 1 / 2, 1, 3 / 2, 2, 3)
  group <- c("Decisive", "Very strong", "Strong", "Substantial", "Barely worth mentioning", "Substantial", "Strong", "Very strong", "Decisive")
  levels <- c("Barely worth mentioning", "Substantial", "Strong", "Very strong", "Decisive")
  highlights <- dplyr::tibble(from = x_from, to = x_to, group = group)
  Zissou1 <- c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")

  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0, color = "darkgray", linetype = "dashed") +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous(
      name = "Exponential                                               log10(BF)                                               Logistic",
      breaks = c(-3, -2, -1, 0, 1, 2, 3),
      labels = c("Inf", "2", "1", "0", "1", "2", "Inf"),
      limits = c(-3, 3)
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(10, "pt"), ends = "both"))
    ) +
    ggplot2::geom_rect(
      data = highlights,
      ggplot2::aes(
        xmin = .data$from,
        xmax = .data$to,
        ymin = -Inf,
        ymax = Inf,
        fill = factor(.data$group, levels = levels)
      ),
      alpha = .5
    ) +
    ggplot2::scale_fill_manual(values = Zissou1) +
    ggplot2::labs(fill = "") +
    ggplot2::theme(legend.position = "bottom")
}

empty_bayes_factor_plot <- function() {
  max_value <- 10

  x_from <- seq(-max_value, max_value - 1, by = 1)
  x_from <- x_from[x_from != 0]
  x_to <- seq(-max_value + 1, max_value, by = 1)
  x_to <- x_to[x_to != 0]
  group <- abs(seq(-max_value + 1, max_value - 1, by = 1))
  highlights <- dplyr::tibble(from = x_from, to = x_to, group = group)

  palette <- c("white", "#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#0868ac", "#084081")

  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_rect(
      data = highlights,
      ggplot2::aes(
        xmin = .data$from,
        xmax = .data$to,
        ymin = -Inf,
        ymax = Inf,
        fill = factor(.data$group)
      ),
      alpha = .5
    ) +
    ggplot2::geom_vline(xintercept = 0, color = "darkgray", linetype = "dashed") +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous(
      name = "log10(BF)",
      breaks = seq(-max_value, max_value, by = 1),
      labels = abs(seq(-max_value, max_value, by = 1)),
      limits = c(-max_value, max_value)
    ) +
    ggplot2::scale_y_continuous(
      name = "Exponential",
      sec.axis = ggplot2::sec_axis(trans = ~ . * 1, name = "Logistic")
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(10, "pt"), ends = "both"))
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(fill = "", y = "pippo") +
    ggplot2::theme(legend.position = "none")
}


#' Plot Omega Values from Mixture Model Fit
#'
#' Generates a plot of the omega values obtained from fitting a mixture model to the data in a `bipod` object.
#' The plot type can be customized to show a histogram, violin plot, or boxplot.
#'
#' @param x A `bipod` object. Must contain a 'fit' field and the 'omega_mixture_model' field in its metadata.
#' @param plot_type A character string specifying the type of plot to generate.
#'  Options are "hist" for histogram, "violin" for violin plot, and "boxplot" for boxplot. (default is 'hist')
#' @param color A character string specifying the color to use for the plot elements. (default is 'maroon')
#'
#' @return A `ggplot2` object displaying the plot of omega values.
#'

#' @export
plot_mixture_model_omega <- function(x, plot_type = "hist", color = "maroon") {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")
  if (!("omega_mixture_model" %in% names(x$metadata))) stop("Input must contain a 'omega_mixture_model' field in its metadata!")
  if (!(plot_type %in% c("hist", "violin", "boxplot"))) stop('plot_type must be one of "hist", "violin", "boxplot"')

  omega_samples <- x$metadata$omega_mixture_model %>%
    dplyr::as_tibble() %>%
    tidyr::pivot_longer(cols = c("omega"))

  if (plot_type == "hist") {
    p <- ggplot2::ggplot(data = omega_samples, mapping = ggplot2::aes(x = .data$value, y = ggplot2::after_stat(..density..))) +
      ggplot2::geom_histogram(fill = color, col = color, alpha = .6, bins = 100) +
      ggplot2::geom_density(col = color, linewidth = 1.2) +
      my_ggplot_theme() +
      ggplot2::scale_x_continuous(
        name = expression(omega),
        breaks = c(0, 1),
        labels = c("Logistic", "Exponential"),
        limits = c(0, 1)
      ) +
      ggplot2::geom_vline(xintercept = 0.5, color = "darkgray", linetype = "dashed")
  }

  if (plot_type == "violin") {
    p <- ggplot2::ggplot(data = omega_samples, mapping = ggplot2::aes(x = .data$name, y = .data$value, fill = .data$name)) +
      ggplot2::geom_violin() +
      ggplot2::lims(y = c(0, 1)) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c(ggplot2::alpha(color, .8))) +
      my_ggplot_theme() +
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
    p <- ggplot2::ggplot(data = omega_samples, mapping = ggplot2::aes(x = .data$name, y = .data$value, fill = .data$name)) +
      ggplot2::geom_boxplot() +
      ggplot2::lims(y = c(0, 1)) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c(ggplot2::alpha(color, .8))) +
      my_ggplot_theme() +
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
