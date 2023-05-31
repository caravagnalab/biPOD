
#' Plot the logarithm of the ELBO and the delta ELBO mean
#' obtained during the variational sampling.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param fit Fit of x from which the elbo information must be extracted
#' @param diagnose A Boolean indicating whether the plots should be colored and
#'  contain info regarding the convergence of the variational sampling.
#'
#' @return A ggplot object with traces of the specified parameters
#' @export
plot_elbo = function(x, fit, diagnose = TRUE) {
  if (!(inherits(x, "bipod"))) stop("The input 'x' must be a 'bipod' object")
  if (!(x$metadata$sampling == "variational")) stop("'plot_elbo' accepts only biPOD objects that have been fitted using variational inference")

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
    qc = line_color = "darkorange2"
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
}

