
#' Plot the fit over the input data.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'

#' @returns A plot of the fit over the input data.
#' @export
plot_fit = function(x) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fits" %in% names(x))) stop("Input must contain a 'fits' field")

  growth_type <- x$fit_info$growth_type

  if (growth_type == "exponential") {
    p <- plot_exponential_fit(x)

  } else {
    p <- plot_logistic_fit(x)
  }
  return(p)
}

plot_exponential_fit = function(x) {

  p <- ggplot2::ggplot()
  data <- data.frame()
  groups <- unique(x$counts$group)

  for (i in 2:length(groups)) {
    previous <- x$counts %>%
      dplyr::filter(.data$group == groups[i-1])

    previous_t <- previous$time[nrow(previous)]
    previous_n <- previous$count[nrow(previous)]

    current <- x$counts %>%
      dplyr::filter(.data$group == groups[i])

    final_t <- current$time[nrow(current)]
    xs <- seq(0, final_t - previous_t, length = 100)

    # Extract fit info
    fit_name <- paste0("fit", groups[i])
    fit <- x$fits[[fit_name]]
    ros <- unname(stats::quantile(rstan::extract(fit, pars="ro")$ro, c(.05, .5, .95)))

    N_low = previous_n * exp(ros[1] * xs)
    N_medium = previous_n * exp(ros[2] * xs)
    N_high = previous_n * exp(ros[3] * xs)

    d <- data.frame(x=xs + previous_t, yl=N_low, ym=N_medium, yh=N_high)
    p <- p +
      ggplot2::geom_line(d, mapping=ggplot2::aes(x=.data$x, y=.data$ym), col="forestgreen") +
      ggplot2::geom_ribbon(d, mapping=ggplot2::aes(x=.data$x, y=.data$ym, ymin=.data$yl, ymax=.data$yh), fill="forestgreen", alpha=.5)
  }

  p <- p +
    ggplot2::geom_point(x$counts, mapping=ggplot2::aes(x=.data$time, y=.data$count)) +
    ggplot2::labs(
      title = paste("Exponential fit", x$sample),
      x = "Time",
      y = "Count"
    ) +
    my_ggplot_theme()
  p
}

plot_logistic_fit = function(x) {

  p <- ggplot2::ggplot()
  data <- data.frame()
  groups <- unique(x$counts$group)

  for (i in 2:length(groups)) {
    previous <- x$counts %>%
      dplyr::filter(.data$group == groups[i-1])

    previous_t <- previous$time[nrow(previous)]
    previous_n <- previous$count[nrow(previous)]

    current <- x$counts %>%
      dplyr::filter(.data$group == groups[i])

    final_t <- current$time[nrow(current)]
    xs <- seq(0, final_t - previous_t, length = 100)

    # Extract fit info
    fit_name <- paste0("fit", groups[i])
    fit <- x$fits[[fit_name]]
    ros <- unname(stats::quantile(rstan::extract(fit, pars="ro")$ro, c(.05, .5, .95)))
    K <- mean(rstan::extract(fit, pars="K")$K)

    N_low = K * previous_n / (previous_n + (K - previous_n) * exp(-xs * ros[1]))
    N_medium = K * previous_n / (previous_n + (K - previous_n) * exp(-xs * ros[2]))
    N_high = K * previous_n / (previous_n + (K - previous_n) * exp(-xs * ros[3]))

    d <- data.frame(x=xs + previous_t, yl=N_low, ym=N_medium, yh=N_high)
    p <- p +
      ggplot2::geom_line(d, mapping=ggplot2::aes(x=.data$x, y=.data$ym), col="forestgreen") +
      ggplot2::geom_ribbon(d, mapping=ggplot2::aes(x=.data$x, y=.data$ym, ymin=.data$yl, ymax=.data$yh), fill="forestgreen", alpha=.5)
  }

  p <- p +
    ggplot2::geom_point(x$counts, mapping=ggplot2::aes(x=.data$time, y=.data$count)) +
    ggplot2::labs(
      title = paste("Exponential fit", x$sample),
      x = "Time",
      y = "Count"
    ) +
    my_ggplot_theme()
  p
}
