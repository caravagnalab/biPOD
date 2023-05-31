
#' Plot the fit over the input data.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit' and must have been fitted with model selection.
#'
#' @returns A plot of the Bayes Factor with its significance.
#' @export
plot_bayes_factor = function(x) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")
  if (!("bayes_factor" %in% names(x$metadata))) stop("Input must contain a 'bayes_factor' field in its metadata!")
  if (!("best_growth" %in% names(x$metadata))) stop("Input must contain a 'best_growth' field in its metadata!")

  logBF <- x$metadata$bayes_factor %>% as.numeric() %>% log10()
  logBF <- min(3, logBF)

  if (x$metadata$best_growth == "Exponential") logBF <- logBF * (-1)

  p <- empty_bayes_factor_plot() +
    geom_point(data=NULL, mapping = ggplot2::aes(x=logBF, y = 0), size = 5)
  p
}

empty_bayes_factor_plot <- function() {
  x_from <- c(-3, -2, -3/2, -1, -1/2, 1/2, 1, 3/2, 2)
  x_to <- c(-2, -3/2, -1, -1/2, 1/2, 1, 3/2, 2, 3)
  group <- c("Decisive", "Very strong",  "Strong", "Substantial", "Barely worth mentioning", "Substantial", "Strong", "Very strong", "Decisive")
  levels = c("Barely worth mentioning", "Substantial", "Strong","Very strong",  "Decisive")
  highlights <- dplyr::tibble(from = x_from, to = x_to, group = group)
  Zissou1 = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")

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
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(type='closed', length = unit(10,'pt'), ends = "both"))
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
