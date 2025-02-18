#' Generate a Comprehensive Model Analysis Report
#'
#' @description
#' Creates a multi-panel visualization report that combines key analytical plots into a single
#' comprehensive display. The report can include:
#' * Model fit visualization with confidence intervals
#' * Breakpoint posterior distributions (if breakpoints were inferred)
#' * Growth rate posterior distributions
#' * Model selection diagnostics (if multiple models were compared)
#'
#' Each component is automatically included based on the available data and analysis
#' in the bipod object, creating a tailored summary of the modeling results.
#'
#' @param x A bipod object that must contain:
#'        * 'fit' field with model fitting results
#'        * metadata fields with model selection information
#'        * Optional: breakpoints analysis results
#'        * Optional: ELBO values for variational inference
#' @param fit_type Character string controlling fit plot complexity:
#'        * "simple": Basic fit visualization
#'        * "complex": Detailed fit visualization with additional elements
#'        Default is "complex"
#' @param breakpoints_color Character string specifying color for breakpoint
#'        posterior visualization. Default is "darkgray"
#' @param shadows_colors Character vector defining colors for time window
#'        highlighting in fit plot. Default is NULL
#' @param t0_posterior_color Character string specifying color for t0 (initial time)
#'        posterior visualization. Default is "darkorange"
#' @param full_process Logical indicating whether to show complete process
#'        details in simple fit plots. Only applies when fit_type = "simple".
#'        Default is FALSE
#'
#' @return A patchwork object containing:
#'   * Main panel: Model fit visualization
#'   * Secondary panels (if applicable):
#'     * Breakpoint posterior distributions
#'     * Growth rate posterior distributions
#'     * Model selection results
#'   All panels are arranged vertically with appropriate size ratios
#'
#' @details
#' The function automatically adapts its output based on the analyses present in
#' the input bipod object:
#' * Breakpoint analysis panels only appear if breakpoints were inferred
#' * Model selection panels only appear if multiple models were compared
#' * Panel heights are automatically adjusted for optimal visualization
#'
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

    # cli::cli_alert_info("Preparing breakpoints posterior diagnostic plot...")
    # if ("breakpoints_elbo" %in% names(x)) {
    #   plots$breakpoints_diagnostic <- biPOD::plot_elbo(x, x$breakpoints_elbo, diagnose = F) +
    #     ggplot2::labs(title = "ELBO - breakpoints inference")
    # } else {
    #   pars <- paste0("changing_times[", 1:length(x$metadata$breakpoints), "]")
    #   plots$breakpoints_diagnostic <- biPOD::plot_traces(x, fit = x$breakpoints_fit, pars = pars)
    #   plots$breakpoints_diagnostic <- ggpubr::annotate_figure(plots$breakpoints_diagnostic, top = ggpubr::text_grob("MCMC traces - breakpoints inference"))
    # }
  }

  # Get rho inference
  cli::cli_alert_info("Preparing growth rates posterior plot...")
  plots$growth_rates <- biPOD::plot_normalized_growth_rate_posteriors(x, colors = shadows_colors) +
    ggplot2::labs(x = bquote(rho), y = "normalized density", title = "Growth rates") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_vline(xintercept = 0, color = "black")

  # Get diagnostic
  # cli::cli_alert_info("Preparing growth rates posterior diagnostic plot...")
  # if ("fit_elbo" %in% names(x)) {
  #   plots$fit_diagnostic <- biPOD::plot_elbo(x, x$fit_elbo, diagnose = F) +
  #     ggplot2::labs(title = "ELBO - fit inference")
  # } else {
  #   plots$fit_diagnostic <- biPOD::plot_traces(x, x$fit)
  #   plots$fit_diagnostic <- ggpubr::annotate_figure(plots$fit_diagnostic, top = ggpubr::text_grob("MCMC traces - fit inference"))
  # }

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
    panel_fits <- patchwork::wrap_plots(plots$fit, plots$breakpoints, plots$growth_rates, plots$model_selection, ncol = 1, heights = c(2, 1, 1, .5))
  } else {
    panel_fits <- patchwork::wrap_plots(plots$fit, plots$growth_rates, plots$model_selection, ncol = 1, heights = c(2, 1, .5))
  }

  panel_fits
}
