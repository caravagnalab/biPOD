#' Plot the fit over the input data.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit' and must have been fitted with model selection.
#' @param with_categories Boolean. If TRUE it will plot the bayes factor according
#' to the categories proposed by Jeffreys in 'The Theory of Probability'
#'
#' @returns A plot of the Bayes Factor with its significance.
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
      ggplot2::geom_point(data = NULL, mapping = ggplot2::aes(x = logBF, y = 0), size = 5)
  } else {
    p <- empty_bayes_factor_plot() +
      ggplot2::geom_point(data = NULL, mapping = ggplot2::aes(x = logBF, y = 0), size = 5)
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
      name = "Exponential                                               log10(BF)                                               Logistic",
      breaks = seq(-max_value, max_value, by = 1),
      labels = abs(seq(-max_value, max_value, by = 1)),
      limits = c(-max_value, max_value)
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(10, "pt"), ends = "both"))
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(fill = "") +
    ggplot2::theme(legend.position = "none")
}
