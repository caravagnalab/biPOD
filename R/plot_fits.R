
#' Plot the fit over the input data.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'

#' @returns A plot of the fit over the input data.
#' @export
fit_plot = function(x) {
  # Check input
  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  assertthat::assert_that("fit" %in% names(x), msg = "Input must contain a 'fit' field")

  growth_type <- x$fit_info$growth_type
  fix_rates <- x$fit_info$fix_rates

  if (growth_type == "exponential") {
    p <- plot_exponential_fit(x, fix_rates)

  } else {
    p <- plot_logistic_fit(x, fix_rates)
  }

  p <- p +
    ggplot2::labs(
      x = "time",
      y = "count"
    )

  return(p)
}

plot_exponential_fit = function(x, fix_rates) {

  d <- x$counts
  n0 <- x$counts$count[1] / x$fit_info$factor_size
  dt <- max(x$counts$time) / 500

  if (fix_rates == 0) {
    n_ros <- 1
    delta_i <- nrow(x$counts) - 1
  } else {
    n_ros <- nrow(x$counts) - 1
    delta_i <- 1
  }

  # Compute ros summary
  ros_df <- purrr::map_df(1:n_ros, ~{
    ros <- as.numeric(stats::quantile(rstan::extract(x$fit, pars=paste0("ro[", .x, "]"))$ro, c(.05, .5, .95)))
    data.frame(low = ros[1], medium = ros[2], high = ros[3])
  })

  # create plot
  p <- ggplot2::ggplot()

  # Compute upper, median and lower bound for predicted counts
  # for every different value of ros
  for (i in 1:n_ros) {
    current_t = d$time[i + delta_i]
    previous_n = d$count[i]

    if (i == 1) {
      times <- seq(0, current_t, by = dt)
    } else {
      times <- seq(0, current_t - d$time[i], by = dt)
    }

    N_low = previous_n * exp(ros_df$low[i] * times)
    N_medium = previous_n * exp(ros_df$medium[i] * times)
    N_high = previous_n * exp(ros_df$high[i] * times)

    p <- p +
      ggplot2::geom_line(data.frame(x=times + d$time[i], y=N_medium), mapping=ggplot2::aes(x=.data$x, y=.data$y), col="forestgreen") +
      ggplot2::geom_ribbon(data.frame(x=times + d$time[i], y=N_medium, yl=N_low, yh=N_high), mapping=ggplot2::aes(x=.data$x, y=.data$y, ymin=.data$yl, ymax=.data$yh), fill="forestgreen", alpha=.5)

  }

  p <- p +
    ggplot2::geom_point(d, mapping=ggplot2::aes(x=.data$time, y=.data$count)) +
    ggplot2::ggtitle(paste("Exponential fit", x$sample)) +
    my_ggplot_theme()

  return(p)
}

plot_logistic_fit = function(x, fix_rates) {

  d <- x$counts
  n0 <- x$counts$count[1] / x$fit_info$factor_size
  K <- mean(rstan::extract(x$fit, pars="K")$K)
  dt <- max(x$counts$time) / 500

  if (fix_rates == 0) {
    n_ros <- 1
    delta_i <- nrow(x$counts) - 1
  } else {
    n_ros <- nrow(x$counts) - 1
    delta_i <- 1
  }

  # Compute ros summary
  ros_df <- purrr::map_df(1:n_ros, ~{
    ros <- as.numeric(stats::quantile(rstan::extract(x$fit, pars=paste0("ro[", .x, "]"))$ro, c(.05, .5, .95)))
    data.frame(low = ros[1], medium = ros[2], high = ros[3])
  })

  # create plot
  p <- ggplot2::ggplot()

  # Compute upper, median and lower bound for predicted counts
  # for every different value of ros
  for (i in 1:n_ros) {
    current_t = d$time[i + delta_i]
    previous_n = d$count[i]

    if (i == 1) {
      times <- seq(0, current_t, by = dt)
    } else {
      times <- seq(0, current_t - d$time[i], by = dt)
    }

    N_low = K * previous_n / (previous_n + (K - previous_n) * exp(-times*(ros_df$low[i])))
    N_medium = K * previous_n / (previous_n + (K - previous_n) * exp(-times*(ros_df$medium[i])))
    N_high = K * previous_n / (previous_n + (K - previous_n) * exp(-times*(ros_df$high[i])))

    p <- p +
      ggplot2::geom_line(data.frame(x=times + d$time[i], y=N_medium), mapping=ggplot2::aes(x=.data$x, y=.data$y), col="forestgreen") +
      ggplot2::geom_ribbon(data.frame(x=times + d$time[i], y=N_medium, yl=N_low, yh=N_high), mapping=ggplot2::aes(x=.data$x, y=.data$y, ymin=.data$yl, ymax=.data$yh), fill="forestgreen", alpha=.5)
  }

  p <- p +
    ggplot2::geom_point(d, mapping=ggplot2::aes(x=.data$time, y=.data$count)) +
    ggplot2::ggtitle(paste("Exponential fit", x$sample)) +
    my_ggplot_theme()

  return(p)
}
