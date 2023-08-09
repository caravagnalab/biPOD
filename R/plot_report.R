#' Plot traces of specified parameters from a fitted Stan model
#'
#' @param x A biPOD object of class `bipod`. Must contains at least one fit
#' @param fit_type .
#' @param breakpoints_color .
#' @param shadows_colors .
#' @param t0_posterior_color .
#' @param full_process .
#'
#' @return A ggplot object with traces of the specified parameters
#' @export
plot_report <- function(x,
                        fit_type = "complex",
                        breakpoints_color = "darkgray",
                        shadows_colors = NULL,
                        t0_posterior_color = "darkorange",
                        full_process = FALSE) {
  plots <- list()

  cli::cli_alert_info("Creating report...")

  # Plot of the fit
  cli::cli_alert_info("Preparing fit plot...")
  if (fit_type != "simple") {
    plots$fit <- plot_fit(x,
      shadows_colors = shadows_colors,
      t0_posterior_color = t0_posterior_color
    ) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(title = "Fit", x = "time (year)")
  } else {
    plots$fit <- plot_fit(x,
      full_process = full_process,
      zoom = FALSE,
      shadows_colors = shadows_colors,
      t0_posterior_color = t0_posterior_color
    ) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(title = "Fit", x = "time (year)")
  }

  # Get breakpoints inference
  if ("breakpoints_fit" %in% names(x)) {
    cli::cli_alert_info("Preparing breakpoints posterior plot...")
    plots$breakpoints <- biPOD::plot_breakpoints_posterior(x, colors = rep(breakpoints_color, length(x$metadata$breakpoints))) +
      ggplot2::lims(x = c(min(x$counts$time), max(x$counts$time))) +
      ggplot2::labs(x = "time (year)", y = "density", title = "Inferred breakpoints") +
      ggplot2::theme(legend.position = "none")

    cli::cli_alert_info("Preparing breakpoints posterior diagnostic plot...")
    if ("breakpoints_elbo" %in% names(x)) {
      plots$breakpoints_diagnostic <- biPOD::plot_elbo(x, x$breakpoints_elbo, diagnose = F) +
        ggplot2::labs(title = "ELBO - breakpoints inference")
    } else {
      pars <- paste0("changing_times[", 1:length(x$metadata$breakpoints), "]")
      plots$breakpoints_diagnostic <- biPOD::plot_traces(x, fit = x$breakpoints_fit, pars = pars)
      plots$breakpoints_diagnostic <- ggpubr::annotate_figure(plots$breakpoints_diagnostic, top = ggpubr::text_grob("MCMC traces - breakpoints inference"))
    }
  }

  # Get rho inference
  cli::cli_alert_info("Preparing growth rates posterior plot...")
  plots$growth_rates <- biPOD::plot_normalized_growth_rate_posteriors(x, colors = shadows_colors) +
    ggplot2::labs(x = bquote(rho), y = "normalized density", title = "Growth rates") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_vline(xintercept = 0, color = "black")

  # Get diagnostic
  cli::cli_alert_info("Preparing growth rates posterior diagnostic plot...")
  if ("fit_elbo" %in% names(x)) {
    plots$fit_diagnostic <- biPOD::plot_elbo(x, x$fit_elbo, diagnose = F) +
      ggplot2::labs(title = "ELBO - fit inference")
  } else {
    plots$fit_diagnostic <- biPOD::plot_traces(x, x$fit)
    plots$fit_diagnostic <- ggpubr::annotate_figure(plots$fit_diagnostic, top = ggpubr::text_grob("MCMC traces - fit inference"))
  }

  # Get model selection plot
  if (!(is.null(x$metadata$model_selection_algo))) {
    cli::cli_alert_info("Preparing model selection plot...")
    if (x$metadata$model_selection_algo == "bayes_factor") {
      plots$model_selection <- biPOD::plot_bayes_factor(x) +
        ggplot2::labs(title = "Model selection")
    } else {
      plots$model_selection <- biPOD::plot_mixture_model_omega(x, plot_type = "hist", color = "maroon") +
        ggplot2::labs(title = "Model selection")
    }
  }

  if ("breakpoints_fit" %in% names(x)) {
    panel_fits <- plots$fit / plots$breakpoints / plots$growth_rates / plots$model_selection +
      patchwork::plot_layout(heights = c(2, 1, 1, .5))
    panel_diagnostics <- plots$breakpoints_diagnostic / plots$fit_diagnostic
  } else {
    panel_fits <- plots$fit / plots$growth_rates / plots$model_selection +
      patchwork::plot_layout(heights = c(2, 1, .5))
    panel_diagnostics <- plots$fit_diagnostic
  }

  list(fits = panel_fits, diagnostics = panel_diagnostics)
}
