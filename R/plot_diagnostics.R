
#' Plot traces of specified parameters from a fitted Stan model
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param pars A character vector of parameters to plot.
#' @param regex_pars A character vector of regular expressions to match against the parameter names in the 'fit' object. At least one of 'pars' or 'regex_pars' should contain input.
#' @param n_col Number of columns of plots in the final output. Default is 4.
#'
#' @return A ggplot object with traces of the specified parameters
#' @export
plot_traces = function(x, pars = c(), regex_pars = c(), n_col = 4) {
  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  assertthat::assert_that("fit" %in% names(x), msg = "Input must contain a 'fit' field")
  assertthat::assert_that(length(pars) > 0 | length(regex_pars) > 0, msg = "Either 'pars' or 'regex_pars' should contain at least one input")

  par_list <- c()
  if (length(regex_pars) > 0) {
    for (p in regex_pars) {
      par_list <- c(par_list, names(x$fit)[grep(p, names(x$fit))])
    }
  } else {
    par_list <- c(c(pars))
  }
  par_list <- unique(par_list)

  assertthat::assert_that(length(par_list) > 0, msg = "No parameters found with given input")
  for (p in par_list) {
    assertthat::assert_that(p %in% names(x$fit), msg = "Parameters passed as input are not present in the 'fit' object")
  }

  plots <- lapply(par_list, function(par) {
    rhat <- rstan::Rhat(as.array(x$fit)[,,par])
    chains <- as.data.frame(rstan::extract(x$fit, pars = par, permuted = FALSE))
    names(chains) <- paste0("chain", 1:ncol(chains))
    chains$index <- 1:nrow(chains)
    chains <- reshape2::melt(chains, id.vars = c("index")) %>%
      dplyr::mutate(parameter = par)

    rhat_color <- if(rhat <= 1.01) "forestgreen" else "indianred"

    p <- ggplot2::ggplot(chains, ggplot2::aes(x=.data$index, y=.data$value, col=.data$variable)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ .data$parameter) +
      ggplot2::scale_color_brewer(palette = "Greens", direction = -1) +
      my_ggplot_theme() +
      ggplot2::theme(legend.position = "none",
                     strip.background = ggplot2::element_rect(colour=rhat_color,
                                                              fill=rhat_color))
  })

  n_col <- if(length(par_list) > n_col) n_col else length(par_list)
  plots <- ggpubr::ggarrange(plotlist = plots,
                             ncol = n_col)
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
