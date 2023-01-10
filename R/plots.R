
COLORS <- c("darkgreen", "dodgerblue", "darkorange", "darkorchid3")
base_color <- "#008080"

#' Plot the evolution of the population over time
#'
#' @param x A biPOD object of class `bipod`.
#' @returns A plot. Represents the evolution of the population over time.
#' @import ggplot2
#' @export
evolution_plot <- function(x) {
  # Check input
  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")

  data <- data.frame(x = x$counts$time, y = x$counts$count)

  pop.plot <- ggplot(data, aes(x=.data$x, y=.data$y)) +
    geom_line(col=base_color) +
    geom_point(col=base_color) +
    ggtitle(x$sample) +
    my_ggplot_theme()

  pop.plot
}

#' Plot the posterior distribution for a specific group of parameters.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param p_name A string indicating the desired parameters.
#'
#' @returns A plot. Represents the posterior distribution for the desired parameters.
#' @import ggplot2
#' @importFrom bayesplot mcmc_areas_ridges
#' @export
posterior_plot <- function(x, p_name) {
  # Check input
  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  assertthat::assert_that(p_name %in% c("lambda", "mu", "ro"), msg = "p_name must be one of: lambda, mu, ro")
  assertthat::assert_that("fit" %in% names(x), msg = "Input must contain a 'fit' field")

  if (p_name %in% c("lambda")) {
    title = "Posterior density for birth rates"
  } else if (p_name %in% c("mu")) {
    title = "Posterior density for death rates"
  } else {
    title = "Posterior density for growth rates"
  }

  p <- mcmc_areas_ridges(as.matrix(x$fit), regex_pars = p_name, prob = 1) +
    ggtitle(title) +
    xlab("Value") +
    my_ggplot_theme()
  p
}

#' Plot the posterior predictive checks for counts data.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param ptype A string representing the type of visualization.
#'
#' `intervals` vertical bars with points indicating generated values medians and darker points indicating observed values
#' `ribbon` ribbon of connected intervals with a line through the median of generated values and a darker line connecting observed values
#'
#' @param prob,prob_outer Values between 0 and 1 indicating the desired probability mass to include in the inner and outer intervals. The defaults are prob=0.5 and prob_outer=0.9.
#'
#' @returns A plot. Represents the posterior predictive checks.
#' @import ggplot2
#' @importFrom rstan extract
#' @importFrom bayesplot ppc_intervals ppc_ribbon
#' @export
ppc_plot = function(x, ptype, prob = 0.5, prob_outer = 0.9) {
  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  assertthat::assert_that("fit" %in% names(x), msg = "Input must contain a 'fit' field")
  assertthat::assert_that(ptype %in% c("intervals", "ribbon"), msg = "ptype must be one of: intervals, ribbon")
  assertthat::assert_that(prob >= 0 && prob < 1, msg = "prob should be between 0 and 1 (with 1 excluded)")
  assertthat::assert_that(prob_outer >= 0 && prob_outer <= 1, msg = "prob_outer should be between 0 and 1")
  assertthat::assert_that(length(grep("N_rep", names(x$fit))) >= 1, msg = "ppc_plot not support ppc_plot does not support the type of model used")

  y <- x$counts$count[2:length(x$counts$count)]
  yrep <- extract(x$fit, pars="N_rep")$N_rep

  if (ptype == "intervals") {
    p <- ppc_intervals(y, yrep, prob = prob, prob_outer = prob_outer)
  } else {
    p <- ppc_ribbon(y, yrep, prob = prob, prob_outer = prob_outer, y_draw = "points")
  }

  p <- p +
    xlab("Time") +
    ylab("Count") +
    my_ggplot_theme()
  p
}


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

  } else {
    # Use different growth rates
    S <- nrow(x$counts) - 1
    ros_df <- purrr::map_df(1:S, ~{
      ros <- as.numeric(quantile(extract(x$fit, pars=paste0("ro[", .x, "]"))$ro, c(.05, .5, .95)))
      data_frame(low = ros[1], medium = ros[2], high = ros[3])
    })

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

      N_low = c(N_low, previous_n * exp(ros_df$low[i] * times))
      N_medium = c(N_medium, previous_n * exp(ros_df$medium[i] * times))
      N_high = c(N_high, previous_n * exp(ros_df$high[i] * times))

      Ts = c(Ts, times + d$time[i])
    }

    times + d$time[i]
  }

  # Create plot
  p <- ggplot() +
    geom_point(d, mapping=aes(x=.data$time, y=.data$count)) +
    geom_line(data.frame(x=Ts, y=N_medium), mapping=aes(x=.data$x, y=.data$y), col=base_color) +
    geom_ribbon(data.frame(x=Ts, y=N_medium, yl=N_low, yh=N_high), mapping=aes(x=.data$x, y=.data$y, ymin=.data$yl, ymax=.data$yh), fill="#008080", alpha=.5) +
    ggtitle(paste("Exponential fit", x$sample)) +
    my_ggplot_theme()
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
      data_frame(low = ros[1], medium = ros[2], high = ros[3])
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
    geom_line(data.frame(x=Ts, y=N_medium), mapping=aes(x=.data$x, y=.data$y), col=base_color) +
    geom_ribbon(data.frame(x=Ts, y=N_medium, yl=N_low, yh=N_high), mapping=aes(x=.data$x, y=.data$y, ymin=.data$yl, ymax=.data$yh), fill="#008080", alpha=.5) +
    ggtitle(paste("Logistic fit", x$sample)) +
    my_ggplot_theme()
  return(p)
}


my_ggplot_theme = function() {
  theme_light()
}
