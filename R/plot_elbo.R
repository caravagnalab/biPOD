#' Plot the logarithm of the ELBO and the delta ELBO mean
#' obtained during the variational sampling.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param elbo_data .
#' @param diagnose A Boolean indicating whether the plots should be colored and
#'  contain info regarding the convergence of the variational sampling.
#'
#' @return A ggplot object with traces of the specified parameters
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
