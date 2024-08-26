#' Plot Traces of Specified Parameters from a Fitted Stan Model
#'
#' Generates trace plots for specified parameters from a fitted Stan model contained within a `bipod` object. The plots visualize the sampled values across iterations for each chain.
#'
#' @param x A `bipod` object that contains at least one fitted model.
#' @param fit The fit within the `bipod` object for which the traces should be plotted.
#' @param pars A character vector specifying the names of the parameters to plot. If empty, plots all available parameters.
#' @param diagnose A logical value indicating whether to color the plots based on convergence diagnostics.
#'  If `TRUE`, the plots will include color coding based on the R-hat convergence statistic. (default is FALSE)
#'
#' @return A `ggplot2` object displaying the trace plots for the specified parameters.
#' @export
plot_traces <- function(x, fit, pars = c(), diagnose = FALSE) {
  if (!(inherits(x, "bipod"))) stop("Input 'x' must be a 'bipod' object")
  # if (!(x$metadata$sampling == "mcmc")) stop("'plot_traces' accepts only biPOD objects that have been fitted using MCMC")

  all_pars <- x$fit$parameters
  if (!(length(pars) > 0)) {
    cli::cli_alert_info("The input vector 'pars' is empty. All the following parameters will be reported: {.val {all_pars}}.  It might take some time...", wrap = T)
    pars <- all_pars
  }

  n_chains <- (ncol(fit$draws) / length(all_pars)) %>% as.integer()

  plots <- lapply(pars, function(par) {
    if (!(par %in% all_pars)) {
      cli::cli_abort(c("The available parameters are:  {.val {all_pars}}", "x" = "You passed {.val {par}}"))
    }
    qc <- "forestgreen"

    rhat <- fit$rhat[[par]]

    chains <- fit$draws[,grepl(par, colnames(fit$draws), fixed = T)] %>%
      `colnames<-`(paste0("chain_", 1:n_chains)) %>%
      dplyr::mutate(index = dplyr::row_number()) %>%
      tidyr::pivot_longer(!.data$index, values_to = "value", names_to = "chain") %>%
      dplyr::mutate(parameter = par)

    if (rhat > 1.01) qc <- "indianred"
    if (!diagnose) qc <- "gray"

    p <- ggplot2::ggplot(chains, ggplot2::aes(x = .data$index, y = .data$value, col = .data$chain)) +
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
          ggplot2::labs(subtitle = "Rhat higher than 1.01") +
          ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0))
      } else {
        p <- p +
          ggplot2::labs(subtitle = "Rhat converged") +
          ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0))
      }
    }

    p
  })

  plots <- ggpubr::ggarrange(plotlist = plots)
  plots
}

#' Plot the Logarithm of the ELBO and Delta ELBO Mean
#'
#' Plots the logarithm of the Evidence Lower Bound (ELBO) and the mean change in ELBO (delta ELBO) obtained during variational inference. This helps to assess the convergence and quality of the variational fit.
#'
#' @param x A `bipod` object of class `bipod` that contains the fitted model and ELBO data.
#' @param elbo_data A data frame or tibble containing ELBO values over iterations. Must include columns `iter`, `ELBO`, and `delta_ELBO_mean`.
#' @param diagnose A logical value indicating whether to color the plot based on the convergence of the ELBO.
#'  If `TRUE`, the plot will include color coding based on whether the ELBO converged.
#'
#' @return A `ggplot2` object displaying the ELBO and delta ELBO mean traces.
#' @export
plot_elbo <- function(x, elbo_data, diagnose = TRUE) {
  if (!(inherits(x, "bipod"))) stop("The input 'x' must be a 'bipod' object")

  elbo_converged <- all(elbo_data$convergence)

  if (elbo_converged) {
    qc <- line_color <- "forestgreen"
    msg <- "ELBO converged."
  } else {
    qc <- line_color <- "indianred"
    msg <- "ELBO did not converge."
  }

  if (!diagnose) {
    qc <- "gray"
    line_color <- "black"
  }

  p <- elbo_data %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(ELBO = -log(-.data$ELBO)) %>%
    dplyr::select(.data$iter, .data$ELBO, .data$delta_ELBO_mean) %>%
    reshape2::melt(id.vars = c("iter")) %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(mapping = ggplot2::aes(x = .data$iter, y = -.data$value), col = line_color) +
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
}
