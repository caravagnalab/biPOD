
#' Plot traces of specified parameters from a fitted Stan model
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param pars A character vector of parameters to plot.
#' @param diagnose A Boolean indicating whether the plots should be colored and
#'  contain info regarding the convergence of the MCMC sampling.
#'
#' @return A ggplot object with traces of the specified parameters
#' @export
plot_traces = function(x, pars = c(), diagnose = FALSE) {
  if (!(inherits(x, "bipod"))) stop("Input 'x' must be a 'bipod' object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fit' field. It appears no model has been fitted yet.")
  if (!("fit_info" %in% names(x))) stop("Input must contain a 'fit_info' field")
  if (!(x$fit_info$sampling == "mcmc")) stop("'plot_traces' accepts only biPOD objects that have been fitted using MCMC")

  all_pars <- colnames(as.matrix(x$fit))
  if (!(length(pars) > 0)) {
    cli::cli_alert_info("The input vector 'pars' is empty. All the following parameters will be reported: {.val {all_pars}}.  It might take some time...", wrap = T)
    pars <-all_pars
  }

  plots <- lapply(pars, function(par) {
    if (!(par %in% all_pars)) {
      cli::cli_abort(c("The available parameters are:  {.val {all_pars}}", "x" = "You passed {.val {par}}"))
    }
    qc = "forestgreen"

    rhat <- rstan::Rhat(as.array(x$fit)[,,par])

    chains <- rstan::extract(x$fit, pars = par, permuted = FALSE) %>%
      dplyr::as_tibble() %>%
      `colnames<-`(paste0("chain", 1:ncol(as.array(x$fit)))) %>%
      dplyr::mutate(index = dplyr::row_number()) %>%
      tidyr::pivot_longer(!.data$index, values_to = "value", names_to = "chain") %>%
      dplyr::mutate(parameter = par)

    if (rhat > 1.01) qc = "indianred"
      if (!diagnose) qc = "gray"

    p <- ggplot2::ggplot(chains, ggplot2::aes(x=.data$index, y=.data$value, col=.data$chain)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ .data$parameter, labeller = ggplot2::label_parsed) +
      ggplot2::scale_color_brewer(palette = "Greens", direction = -1) +
      my_ggplot_theme() +
      ggplot2::theme(
        legend.position = "none",
        strip.background = ggplot2::element_rect(fill = qc)
      )

    if (diagnose) {
      if (qc == "indianred") {
        p <- p +
          ggplot2::labs(subtitle ="Rhat higher than 1.01") +
          ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0))
      } else {
        p <- p +
          ggplot2::labs(subtitle ="Rhat converged") +
          ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0))
      }
    }

    p
  })

  plots <- ggpubr::ggarrange(plotlist = plots)
  return(plots)
}

#' Plot the logarithm of the ELBO and the delta ELBO mean
#' obtained during the variational sampling.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param diagnose A Boolean indicating whether the plots should be colored and
#'  contain info regarding the convergence of the variational sampling.
#'
#' @return A ggplot object with traces of the specified parameters
#' @export
plot_elbo = function(x, diagnose = TRUE) {
  if (!(inherits(x, "bipod"))) stop("The input 'x' must be a 'bipod' object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fit' field. It appears no model has been fitted yet.")
  if (!("fit_info" %in% names(x))) stop("Input must contain a 'fit_info' field")
  if (!(x$fit_info$sampling == "variational")) stop("'plot_elbo' accepts only biPOD objects that have been fitted using variational inference")

  elbo_data <- x$elbo_data

  elbo_converged <- all(elbo_data$convergence)
  pareto_k <- elbo_data$pareto_k[1]

  if ((pareto_k > 1) | (!elbo_converged)) {
    qc = line_color = "indianred"
      msg = "Pareto k higher than 1 and/or ELBO not convergent."
  } else if (pareto_k == 0) {
    qc = line_color = "forestgreen"
      msg = "Pareto k lower than 0.5 and convergent ELBO."
  } else {
    qc = line_color = "darkorange"
      msg = "Convergent ELBO but Pareto k between .5 and 1."
  }

  if (!diagnose) {
    qc = "gray"
      line_color = "black"
  }

  p <- elbo_data %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(ELBO = -log(-.data$ELBO)) %>%
    dplyr::select(.data$iter, .data$ELBO, .data$delta_ELBO_mean) %>%
    reshape2::melt(id.vars = c("iter")) %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(mapping = ggplot2::aes(x=.data$iter, y=-.data$value), col = line_color) +
    ggplot2::facet_wrap(~ .data$variable, scales = "free") +
    ggplot2::labs(
      x = "iteration",
      y = "value",
      title = "ELBO"
    ) +
    my_ggplot_theme() +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_rect(fill = qc)
    )

  if (diagnose) {
    p <- p +
      ggplot2::labs(subtitle = msg) +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0))
  }

  return(p)

  #######

  # plots <- lapply(names(x$elbo_data), function(elbo_name) {
  #   elbo_data <- x$elbo_data[[elbo_name]]
  #
  #   elbo_converged <- all(elbo_data$convergence)
  #   pareto_k <- elbo_data$pareto_k[1]
  #
  #   if ((pareto_k > 1) | (!elbo_converged)) {
  #     qc = line_color = "indianred"
  #     msg = "Pareto k higher than 1 and/or ELBO not convergent."
  #   } else if (pareto_k == 0) {
  #     qc = line_color = "forestgreen"
  #     msg = "Pareto k lower than 0.5 and convergent ELBO."
  #   } else {
  #     qc = line_color = "darkorange"
  #     msg = "Convergent ELBO but Pareto k between .5 and 1."
  #   }
  #
  #   if (!diagnose) {
  #     qc = "gray"
  #     line_color = "black"
  #   }
  #
  #   if (add_title) {
  #     title = paste0("Group ", gsub("elbo", "", elbo_name))
  #   } else {
  #     title = ""
  #   }
  #
  #   p <- elbo_data %>%
  #     as.data.frame %>%
  #     dplyr::mutate(ELBO = -log(-.data$ELBO)) %>%
  #     dplyr::select(.data$iter, .data$ELBO, .data$delta_ELBO_mean) %>%
  #     reshape2::melt(id.vars = c("iter")) %>%
  #     ggplot2::ggplot() +
  #     ggplot2::geom_line(mapping = ggplot2::aes(x=.data$iter, y=-.data$value), col = line_color) +
  #     ggplot2::facet_wrap(~ .data$variable, scales = "free") +
  #     ggplot2::labs(
  #       x = "iteration",
  #       y = "value",
  #       title = title
  #     ) +
  #     my_ggplot_theme() +
  #     ggplot2::theme(
  #       legend.position = "none",
  #       strip.background = ggplot2::element_rect(fill = qc)
  #     )
  #
  #   if (diagnose) {
  #     p <- p +
  #       ggplot2::labs(subtitle = msg) +
  #       ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0))
  #   }
  #
  #   p
  # })
  # plots <- ggpubr::ggarrange(plotlist = plots, ncol=1)
  # plots
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
ppc_plot = function(x, ptype, prob = 0.5, prob_outer = 0.9) {
  cli::cli_abort("TODO")
  # assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  # assertthat::assert_that("fit" %in% names(x), msg = "Input must contain a 'fit' field")
  # assertthat::assert_that(ptype %in% c("intervals", "ribbon"), msg = "ptype must be one of: intervals, ribbon")
  # assertthat::assert_that(prob >= 0 && prob < 1, msg = "prob should be between 0 and 1 (with 1 excluded)")
  # assertthat::assert_that(prob_outer >= 0 && prob_outer <= 1, msg = "prob_outer should be between 0 and 1")
  # assertthat::assert_that(length(grep("N_rep", names(x$fit))) >= 1, msg = "ppc_plot not support ppc_plot does not support the type of model used")
  #
  # y <- x$counts$count[2:length(x$counts$count)]
  # yrep <- rstan::extract(x$fit, pars="N_rep")$N_rep
  #
  # if (ptype == "intervals") {
  #   p <- bayesplot::ppc_intervals(y, yrep, prob = prob, prob_outer = prob_outer)
  # } else {
  #   p <- bayesplot::ppc_ribbon(y, yrep, prob = prob, prob_outer = prob_outer, y_draw = "points")
  # }
  #
  # p <- p +
  #   ggplot2::labs(
  #     y = 'time',
  #     x = "count"
  #   ) +
  #   my_ggplot_theme()
  # p
}
