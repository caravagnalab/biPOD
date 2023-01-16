
#' Plot traces of specified parameters from a fitted Stan model
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param pars A character vector of parameters to plot.
#'
#' @return A ggplot object with traces of the specified parameters
#' @export
plot_traces = function(x, pars = c()) {
  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  assertthat::assert_that("fits" %in% names(x), msg = "Input must contain a 'fits' field")
  assertthat::assert_that(length(pars) > 0, msg = "'pars' should contain at least one input")

  plots <- lapply(names(x$fits), function(fit_name) {
    fit <- x$fits[[fit_name]]
    long_chains <- data.frame()
    for (par in pars) {
      rhat <- rstan::Rhat(as.array(fit)[,,par])
      chains <- as.data.frame(rstan::extract(fit, pars = par, permuted = FALSE))
      names(chains) <- paste0("chain", 1:ncol(chains))
      chains$index <- 1:nrow(chains)
      chains <- reshape2::melt(chains, id.vars = c("index")) %>%
        dplyr::mutate(parameter = paste0(par, gsub("fit", "", fit_name))) %>%
        dplyr::mutate(color = ifelse(rhat <= 1.01, "forestgreen", "indianred"))

      long_chains <- dplyr::bind_rows(long_chains, chains)
    }

    p <- ggplot2::ggplot(long_chains, ggplot2::aes(x=.data$index, y=.data$value, col=.data$variable)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ .data$parameter) +
      ggplot2::scale_color_brewer(palette = "Greens", direction = -1) +
      my_ggplot_theme() +
      ggplot2::theme(legend.position = "none")
    p
  })

  plots <- ggpubr::ggarrange(plotlist = plots, ncol=1)
  return(plots)
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
#' @export
ppc_plot = function(x, ptype, prob = 0.5, prob_outer = 0.9) {
  cli::cli_alert_danger("TODO")
  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  assertthat::assert_that("fit" %in% names(x), msg = "Input must contain a 'fit' field")
  assertthat::assert_that(ptype %in% c("intervals", "ribbon"), msg = "ptype must be one of: intervals, ribbon")
  assertthat::assert_that(prob >= 0 && prob < 1, msg = "prob should be between 0 and 1 (with 1 excluded)")
  assertthat::assert_that(prob_outer >= 0 && prob_outer <= 1, msg = "prob_outer should be between 0 and 1")
  assertthat::assert_that(length(grep("N_rep", names(x$fit))) >= 1, msg = "ppc_plot not support ppc_plot does not support the type of model used")

  y <- x$counts$count[2:length(x$counts$count)]
  yrep <- rstan::extract(x$fit, pars="N_rep")$N_rep

  if (ptype == "intervals") {
    p <- bayesplot::ppc_intervals(y, yrep, prob = prob, prob_outer = prob_outer)
  } else {
    p <- bayesplot::ppc_ribbon(y, yrep, prob = prob, prob_outer = prob_outer, y_draw = "points")
  }

  p <- p +
    ggplot2::labs(
      y = 'time',
      x = "count"
    ) +
    my_ggplot_theme()
  p
}
