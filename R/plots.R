
COLORS <- c("darkgreen", "dodgerblue", "darkorange", "darkorchid3")
base_color <- "#008080"

#' Plot the evolution of the population over time
#'
#' @param x A biPOD object of class `bipod`.
#' @returns A plot. Represents the evolution of the population over time.
#' @import ggplot2
#' @export
evolution_plot <- function(x) {
  stopifnot(inherits(x, 'bipod'))

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
  stopifnot(inherits(x, "bipod"))
  stopifnot(p_name %in% c("lambda", "mu", "ro"))
  stopifnot('fit' %in% names(x))

  if (p_name %in% c("lambda")) {
    title = "Posterior density for birth rates"
  } else if (p_name %in% c("mu")) {
    title = "Posterior density for death rates"
  } else {
    title = "Posterior density for growth rates"
  }

  p <- bayesplot::mcmc_areas_ridges(as.matrix(x$fit), regex_pars = p_name, prob = 1) +
    ggplot2::ggtitle(title) +
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
#' @importFrom bayesplot ppc_intervals ppc_ribbon
#' @export
ppc_plot = function(x, ptype, prob = 0.5, prob_outer = 0.9) {
  stopifnot(ptype %in% c("intervals", "ribbon"))
  stopifnot(prob >= 0 && prob <= 1)
  stopifnot(prob_outer >= 0 && prob <= 1)

  y <- x$counts$count[2:length(x$counts$count)]
  yrep <- rstan::extract(x$fit, pars="N_rep")$N_rep

  if (ptype == "intervals") {
    p <- bayesplot::ppc_intervals(y, yrep, prob = prob, prob_outer = prob_outer)
  } else {
    p <- bayesplot::ppc_ribbon(y, yrep, prob = prob, prob_outer = prob_outer, y_draw = "points")
  }

  p <- p + ggplot2::xlab("Time") + ggplot2::ylab("Count") + my_ggplot_theme()
  p
}


#' Plot the fit over the input data.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'

#' @returns A plot.
#' @import ggplot2
#' @export
fit_plot = function(x) {
  stopifnot(inherits(x, "bipod"))
  stopifnot('fit' %in% names(x))

  growth_type <- x$fit_info$growth_type
  fix_rates <- x$fit_info$fix_rates

  if (growth_type == "exponential") {

    if (fix_rates == 0) {
      ros <- as.numeric(quantile(rstan::extract(x$fit, pars="ro")$ro, c(.05, .5, .95)))

      d <- x$counts
      n0 <- x$counts$count[1] / x$fit_info$factor_size

      Ts <- seq(0, max(d$time), by=0.1)
      N_low <- n0 * exp(ros[1] * Ts)
      N_medium <- n0 * exp(ros[2] * Ts)
      N_high <- n0 * exp(ros[3] * Ts)

      p <- ggplot2::ggplot() +
        ggplot2::geom_point(d, mapping=aes(x=time, y=count)) +
        ggplot2::geom_line(data.frame(x=Ts, y=N_medium), mapping=aes(x=x, y=y), col=base_color) +
        ggplot2::geom_ribbon(data.frame(x=Ts, y=N_medium, yl=N_low, yh=N_high), mapping=aes(x=x, y=y, ymin=yl, ymax=yh), fill="#008080", alpha=.5) +
        ggtitle("Exponential fit") +
        my_ggplot_theme()
    } else {

      # Use all ros
      S <- nrow(x$counts) - 1
      ros_medium = rep(0, S)
      ros_low = ros_high = ros_medium
      for (i in 1:S) {
        ros <- as.numeric(quantile(rstan::extract(x$fit, pars=paste0("ro[", i, "]"))$ro, c(.05, .5, .95)))
        ros_low[i] = ros[1]
        ros_medium[i] = ros[2]
        ros_high[i] = ros[3]
      }

      d <- x$counts
      n0 <- x$counts$count[1] / x$fit_info$factor_size

      Ts = c()
      N_medium = c()
      N_low = N_high = N_medium

      for (i in 1:S) {
        current_t = d$time[i+1]
        previous_n = d$count[i]

        if (i == 1) {
          times <- seq(0, current_t, by = 0.1)
        } else {
          times <- seq(0, current_t - d$time[i], by = 0.1)
        }

        N_low = c(N_low, previous_n * exp(ros_low[i] * times))
        N_medium = c(N_medium, previous_n * exp(ros_medium[i] * times))
        N_high = c(N_high, previous_n * exp(ros_high[i] * times))

        Ts = c(Ts, times + d$time[i])
      }

      times + d$time[i]

      p <- ggplot2::ggplot() +
        ggplot2::geom_point(d, mapping=aes(x=time, y=count)) +
        ggplot2::geom_ribbon(data.frame(x=Ts, y=N_medium, yl=N_low, yh=N_high), mapping=aes(x=x, y=y, ymin=yl, ymax=yh), fill="#008080", alpha=.5) +
        ggplot2::ggtitle("Logistic fit") +
        my_ggplot_theme()

    }
  } else {

    if (fix_rates == 0) {
      ros <- as.numeric(quantile(rstan::extract(x$fit, pars="ro")$ro, c(.05, .5, .95)))
      K <- mean(rstan::extract(x$fit, pars="K")$K)

      d <- x$counts
      n0 <- x$counts$count[1] / x$fit_info$factor_size

      Ts <- seq(0, max(d$time), by=0.1)
      N_low <- K * n0 / (n0 + (K - n0) * exp(-Ts*(ros[1])))
      N_medium <- K * n0 / (n0 + (K - n0) * exp(-Ts*(ros[2])))
      N_high <- K * n0 / (n0 + (K - n0) * exp(-Ts*(ros[3])))

      p <- ggplot2::ggplot() +
        ggplot2::geom_point(d, mapping=aes(x=time, y=count)) +
        ggplot2::geom_line(data.frame(x=Ts, y=N_medium), mapping=aes(x=x, y=y), col=base_color) +
        ggplot2::geom_ribbon(data.frame(x=Ts, y=N_medium, yl=N_low, yh=N_high), mapping=aes(x=x, y=y, ymin=yl, ymax=yh), fill="#008080", alpha=.5) +
        ggtitle("Logistic fit") +
        my_ggplot_theme()
    } else {
      S <- nrow(x$counts) - 1
      ros_medium = rep(0, S)
      ros_low = ros_high = ros_medium
      for (i in 1:S) {
        ros <- as.numeric(quantile(rstan::extract(x$fit, pars=paste0("ro[", i, "]"))$ro, c(.05, .5, .95)))
        ros_low[i] = ros[1]
        ros_medium[i] = ros[2]
        ros_high[i] = ros[3]
      }

      K <- mean(rstan::extract(x$fit, pars="K")$K)

      d <- x$counts
      n0 <- x$counts$count[1] / x$fit_info$factor_size

      Ts = c()
      N_medium = c()
      N_low = N_high = N_medium

      for (i in 1:S) {
        current_t = d$time[i+1]
        previous_n = d$count[i]

        if (i == 1) {
          times <- seq(0, current_t, by = 0.1)
        } else {
          times <- seq(0, current_t - d$time[i], by = 0.1)
        }

        N_low = c(N_low, K * previous_n / (previous_n + (K - previous_n) * exp(-times*(ros_low[i]))))
        N_medium = c(N_medium, K * previous_n / (previous_n + (K - previous_n) * exp(-times*(ros_medium[i]))))
        N_high = c(N_high, K * previous_n / (previous_n + (K - previous_n) * exp(-times*(ros_high[i]))))
        Ts = c(Ts, times + d$time[i])
      }

      times + d$time[i]

      p <- ggplot2::ggplot() +
        ggplot2::geom_point(d, mapping=aes(x=time, y=count)) +
        ggplot2::geom_ribbon(data.frame(x=Ts, y=N_medium, yl=N_low, yh=N_high), mapping=aes(x=x, y=y, ymin=yl, ymax=yh), fill="#008080", alpha=.5) +
        ggplot2::ggtitle("Logistic fit") +
        my_ggplot_theme()
    }
  }

  return(p)
}

my_ggplot_theme = function() {
  theme_light()
}
