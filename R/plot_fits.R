
#' Plot the fit over the input data.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'

#' @returns A plot of the fit over the input data.
#'
#' @import ggplot2
#' @importFrom stats quantile
#' @importFrom rstan extract
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
    xlab("Time") +
    ylab("Count")

  return(p)
}

#' @importFrom purrr map_df
plot_exponential_fit = function(x, fix_rates) {

  if (fix_rates == 0) {
    ros <- as.numeric(quantile(extract(x$fit, pars="ro")$ro, c(.05, .5, .95)))

    d <- x$counts
    n0 <- x$counts$count[1] / x$fit_info$factor_size

    Ts <- seq(0, max(d$time), length = 10 * nrow(d))
    N_low <- n0 * exp(ros[1] * Ts)
    N_medium <- n0 * exp(ros[2] * Ts)
    N_high <- n0 * exp(ros[3] * Ts)

    # Create plot
    p <- ggplot() +
      geom_line(data.frame(x=Ts, y=N_medium), mapping=aes(x=.data$x, y=.data$y), col="forestgreen") +
      geom_ribbon(data.frame(x=Ts, y=N_medium, yl=N_low, yh=N_high), mapping=aes(x=.data$x, y=.data$y, ymin=.data$yl, ymax=.data$yh), fill="forestgreen", alpha=.5) +
      geom_point(d, mapping=aes(x=.data$time, y=.data$count)) +
      ggtitle(paste("Exponential fit", x$sample)) +
      my_ggplot_theme()

  } else {
    # Use different growth rates
    S <- nrow(x$counts) - 1
    ros_df <- purrr::map_df(1:S, ~{
      ros <- as.numeric(quantile(extract(x$fit, pars=paste0("ro[", .x, "]"))$ro, c(.05, .5, .95)))
      data.frame(low = ros[1], medium = ros[2], high = ros[3])
    })

    d <- x$counts
    n0 <- x$counts$count[1] / x$fit_info$factor_size

    # Create plot
    p <- ggplot()

    for (i in 1:S) {
      current_t = d$time[i+1]
      previous_n = d$count[i]

      if (i == 1) {
        times <- seq(0, current_t, length=10)
      } else {
        times <- seq(0, current_t - d$time[i], length=10)
      }

      N_low = previous_n * exp(ros_df$low[i] * times)
      N_medium = previous_n * exp(ros_df$medium[i] * times)
      N_high = previous_n * exp(ros_df$high[i] * times)

      p <- p +
        geom_line(data.frame(x=times + d$time[i], y=N_medium), mapping=aes(x=.data$x, y=.data$y), col="forestgreen") +
        geom_ribbon(data.frame(x=times + d$time[i], y=N_medium, yl=N_low, yh=N_high), mapping=aes(x=.data$x, y=.data$y, ymin=.data$yl, ymax=.data$yh), fill="forestgreen", alpha=.5)

    }

    p <- p + geom_point(d, mapping=aes(x=.data$time, y=.data$count)) +
      ggtitle(paste("Exponential fit", x$sample)) +
      my_ggplot_theme()
  }
  p
}

#' @importFrom purrr map_df
plot_logistic_fit = function(x, fix_rates) {
  if (fix_rates == 0) {
    ros <- as.numeric(quantile(extract(x$fit, pars="ro")$ro, c(.05, .5, .95)))
    K <- mean(extract(x$fit, pars="K")$K)

    d <- x$counts
    n0 <- x$counts$count[1] / x$fit_info$factor_size

    Ts <- seq(0, max(d$time), length = 10 * nrow(d))
    N_low <- K * n0 / (n0 + (K - n0) * exp(-Ts*(ros[1])))
    N_medium <- K * n0 / (n0 + (K - n0) * exp(-Ts*(ros[2])))
    N_high <- K * n0 / (n0 + (K - n0) * exp(-Ts*(ros[3])))

  } else {
    # Use all ros
    S <- nrow(x$counts) - 1
    ros_df <- purrr::map_df(1:S, ~{
      ros <- as.numeric(quantile(extract(x$fit, pars=paste0("ro[", .x, "]"))$ro, c(.05, .5, .95)))
      data.frame(low = ros[1], medium = ros[2], high = ros[3])
    })

    K <- mean(extract(x$fit, pars="K")$K)

    d <- x$counts
    n0 <- x$counts$count[1] / x$fit_info$factor_size

    Ts = c()
    N_medium = c()
    N_low = N_high = N_medium

    for (i in 1:S) {
      current_t = d$time[i+1]
      previous_n = d$count[i]

      if (i == 1) {
        times <- seq(0, current_t, length=10)
      } else {
        times <- seq(0, current_t - d$time[i], length=10)
      }

      N_low = c(N_low, K * previous_n / (previous_n + (K - previous_n) * exp(-times*(ros_df$low[i]))))
      N_medium = c(N_medium, K * previous_n / (previous_n + (K - previous_n) * exp(-times*(ros_df$medium[i]))))
      N_high = c(N_high, K * previous_n / (previous_n + (K - previous_n) * exp(-times*(ros_df$high[i]))))
      Ts = c(Ts, times + d$time[i])
    }

    times + d$time[i]
  }

  p <- ggplot() +
    geom_point(d, mapping=aes(x=.data$time, y=.data$count)) +
    geom_line(data.frame(x=Ts, y=N_medium), mapping=aes(x=.data$x, y=.data$y), col="forestgreen") +
    geom_ribbon(data.frame(x=Ts, y=N_medium, yl=N_low, yh=N_high), mapping=aes(x=.data$x, y=.data$y, ymin=.data$yl, ymax=.data$yh), fill="forestgreen", alpha=.5) +
    ggtitle(paste("Logistic fit", x$sample)) +
    my_ggplot_theme()
  return(p)
}
