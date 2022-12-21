
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

  p <- bayesplot::mcmc_areas_ridges(as.matrix(x$fit), regex_pars = "ro", prob = 1) +
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

my_ggplot_theme = function() {
  theme_light()
}
