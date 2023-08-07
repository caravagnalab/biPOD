#' Plot the fit over the input data.
#'
#' @param x A bipod object of class `bipod`. Must contains 'fit'
#' @param CI confidence interval for the growth rate to plot
#' @param add_posteriors Boolean. If TRUE, posteriors regarding
#' the time of origin of the resistant population and the time of death
#' of the sensitive popultation will be added
#'
#' @returns A plot of the fit over the input data.
#' @export
plot_two_pop_fit <- function(x, CI = .95, add_posteriors = TRUE) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("two_pop_fit" %in% names(x))) stop("Input must contain a 'two_pop_fit' field")

  alpha <- 1 - CI
  fitted_data <- get_data_for_two_pop_plot(x, alpha = alpha)

  # Plot data
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(x$counts, mapping = ggplot2::aes(x = .data$time, y = .data$count)) + # original points
    ggplot2::geom_line(fitted_data %>% dplyr::filter(.data$group == "resistant"), mapping = ggplot2::aes(x = .data$x, y = .data$y), col = "forestgreen") +
    ggplot2::geom_ribbon(fitted_data %>% dplyr::filter(.data$group == "resistant"), mapping = ggplot2::aes(x = .data$x, y = .data$y, ymin = .data$ylow, ymax = .data$yhigh), fill = "forestgreen", alpha = .3) +
    ggplot2::geom_line(fitted_data %>% dplyr::filter(.data$group == "sensible"), mapping = ggplot2::aes(x = .data$x, y = .data$y), col = "indianred") +
    ggplot2::geom_ribbon(fitted_data %>% dplyr::filter(.data$group == "sensible"), mapping = ggplot2::aes(x = .data$x, y = .data$y, ymin = .data$ylow, ymax = .data$yhigh), fill = "indianred", alpha = .3) +
    my_ggplot_theme()

  p <- add_posterior(base_plot = p, param = "t0_r", x = x, x_fit = x$two_pop_fit, color = "forestgreen")
  p <- add_posterior(base_plot = p, param = "t_end", x = x, x_fit = x$two_pop_fit, color = "indianred")

  return(p)
}

# Utils

get_data_for_two_pop_plot <- function(x, alpha) {
  fit <- x$two_pop_fit

  factor_size <- x$metadata$factor_size # factor size

  # Produce ro quantiles
  rho_samples <- get_parameters(fit, par_list = c("rho_s", "rho_r"))

  rho_quantiles <- rho_samples %>%
    dplyr::group_by(.data$parameter) %>%
    dplyr::summarise(low = stats::quantile(.data$value, alpha / 2), mid = stats::quantile(.data$value, .5), high = stats::quantile(.data$value, 1 - alpha / 2))
  rho_s_quantiles <- rho_quantiles %>% dplyr::filter(.data$parameter == "rho_s")
  rho_r_quantiles <- rho_quantiles %>% dplyr::filter(.data$parameter == "rho_r")

  median_t0_r <- get_parameter(fit, "t0_r") %>%
    dplyr::pull(.data$value) %>%
    stats::median()

  min_t <- min(median_t0_r, min(x$counts$time))

  xs <- seq(min_t, max(x$counts$time), length = 1000)

  median_ns <- get_parameter(fit, "ns") %>%
    dplyr::pull(.data$value) %>%
    stats::median()

  func <- two_pops_evo

  ylow <- lapply(xs, two_pops_evo, ns = median_ns, t0_r = median_t0_r, rho_s = -rho_s_quantiles$low, rho_r = rho_r_quantiles$low)
  ymid <- lapply(xs, two_pops_evo, ns = median_ns, t0_r = median_t0_r, rho_s = -rho_s_quantiles$mid, rho_r = rho_r_quantiles$mid)
  yhigh <- lapply(xs, two_pops_evo, ns = median_ns, t0_r = median_t0_r, rho_s = -rho_s_quantiles$high, rho_r = rho_r_quantiles$high)

  ylow_r <- lapply(ylow, function(y) {
    y$r_pop
  }) %>% unlist()
  ylow_s <- lapply(ylow, function(y) {
    y$s_pop
  }) %>% unlist()
  ymid_r <- lapply(ymid, function(y) {
    y$r_pop
  }) %>% unlist()
  ymid_s <- lapply(ymid, function(y) {
    y$s_pop
  }) %>% unlist()
  yhigh_r <- lapply(yhigh, function(y) {
    y$r_pop
  }) %>% unlist()
  yhigh_s <- lapply(yhigh, function(y) {
    y$s_pop
  }) %>% unlist()

  r_data <- dplyr::tibble(x = xs, y = factor_size * ymid_r, ylow = factor_size * ylow_r, yhigh = factor_size * yhigh_r, group = "resistant")
  s_data <- dplyr::tibble(x = xs, y = factor_size * ymid_s, ylow = factor_size * ylow_s, yhigh = factor_size * yhigh_s, group = "sensible")

  fitted_data <- dplyr::bind_rows(r_data, s_data)

  return(fitted_data)
}
