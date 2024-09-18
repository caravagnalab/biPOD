#' Plot the Fit Over the Input Data for Two Populations
#'
#' Creates a plot of the fit over the input data for a model with two populations.
#' The plot can display the fit for each population separately and optionally include posteriors for key parameters.
#'
#' @param x A `bipod` object. Must contain a 'two_pop_fit' field.
#' @param CI Numeric value representing the confidence interval for the growth rate to plot. (default is 0.95)
#' @param f_posteriors .
#' @param t_posteriors .
#' @param r_posteriors .
#' @param split_process Logical value indicating whether to plot the dynamics of the two populations separately. (default is FALSE)
#' @param resistant_color .
#' @param sensitive_color .
#'
#' @return A `ggplot2` object showing the fit over the input data.
#' @export
plot_two_pop_fit <- function(
    x,
    CI = .95,
    f_posteriors = TRUE,
    t_posteriors = TRUE,
    r_posteriors = TRUE,
    split_process = TRUE,
    resistant_color = "steelblue",
    sensitive_color = "indianred") {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("two_pop_fit" %in% names(x))) stop("Input must contain a 'two_pop_fit' field")
  #if (add_posteriors & !(split_process)) stop("Posteriors can be only added if 'split_process' = TRUE")

  alpha <- 1 - CI
  fitted_data <- get_data_for_two_pop_plot(x, alpha = alpha)

  fitted_data <- dplyr::bind_rows(fitted_data,
                                  fitted_data %>%
                                    dplyr::select(x, .data$ylow, .data$yhigh, .data$y) %>%
                                    dplyr::group_by(x) %>%
                                    dplyr::summarise_all(sum) %>%
                                    dplyr::mutate(group = 'total'))


  times <- x$counts$time

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(x$counts, mapping = ggplot2::aes(x = .data$time, y = .data$count)) + # original points
    ggplot2::geom_line(fitted_data %>%
                         dplyr::filter(.data$group == "total", .data$x >= min(times), .data$x <= max(times)), mapping = ggplot2::aes(x = .data$x, y = .data$y), col = "black") +
    my_ggplot_theme()

  if (split_process) {
    p <- p +
      ggplot2::geom_line(fitted_data %>% dplyr::filter(.data$group == "resistant", .data$x <= max(times)), mapping = ggplot2::aes(x = .data$x, y = .data$y), col = resistant_color) +
      ggplot2::geom_ribbon(fitted_data %>% dplyr::filter(.data$group == "resistant", .data$x <= max(times)), mapping = ggplot2::aes(x = .data$x, y = .data$y, ymin = .data$ylow, ymax = .data$yhigh), fill = resistant_color, alpha = .3) +
      ggplot2::geom_line(fitted_data %>% dplyr::filter(.data$group == "sensible", .data$x >= min(times)), mapping = ggplot2::aes(x = .data$x, y = .data$y), col = sensitive_color) +
      ggplot2::geom_ribbon(fitted_data %>% dplyr::filter(.data$group == "sensible", .data$x >= min(times)), mapping = ggplot2::aes(x = .data$x, y = .data$y, ymin = .data$ylow, ymax = .data$yhigh), fill = sensitive_color, alpha = .3)
  }

  plots <- list(p = p)

  if (t_posteriors) {
    t_data <- get_parameters(x$two_pop_fit, par_list = c('t0_r', "t_end")) %>%
      dplyr::mutate(parameter = dplyr::if_else(.data$parameter == "t0_r", "t_r", 't_e'))

    time_limits <- c(min(t_data$value, min(fitted_data$x)), max(t_data$value, max(fitted_data$x)))

    p <- p + ggplot2::lims(x = time_limits)

    t_plot <- t_data%>%
      ggplot2::ggplot(mapping = ggplot2::aes(x=.data$value, fill=.data$parameter)) +
      ggplot2::geom_histogram(alpha = .7, binwidth = 0.005) +
      #scale_color_manual(values = c("t_e"=sensitive_color, "t_r"=resistant_color)) +
      ggplot2::scale_fill_manual(values = c("t_e"=sensitive_color, "t_r"=resistant_color)) +
      my_ggplot_theme() +
      ggplot2::labs(fill = "", col="", x="time") +
      ggplot2::lims(x = time_limits)

    plots = list(p = p, t_plot=t_plot)
  }

  if (r_posteriors) {
    d <- get_parameters(x$two_pop_fit, par_list = c('rho_r', "rho_s")) %>%
      dplyr::mutate(value = ifelse(.data$parameter == "rho_s", -.data$value, .data$value))

    rho_plot <- d %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x=.data$value,fill=.data$parameter)) +
      #geom_density(alpha = .3) +
      ggplot2::geom_histogram(alpha = .7, binwidth = 0.005) +
      #scale_color_manual(values = c("rho_s"=sensitive_color, "rho_r"=resistant_color)) +
      ggplot2::scale_fill_manual(values = c("rho_s"=sensitive_color, "rho_r"=resistant_color)) +
      my_ggplot_theme() +
      ggplot2::labs(fill = "", col="", x="Growth rate") +
      ggplot2::geom_vline(xintercept = 0, linetype='dashed')

    plots$rho_plot = rho_plot
  }

  if (f_posteriors) {
    ns <- get_parameters(x$two_pop_fit, par_list = c('ns[1]')) %>% dplyr::pull(.data$value)
    nr <- get_parameters(x$two_pop_fit, par_list = c('nr[1]')) %>% dplyr::pull(.data$value)
    d <- dplyr::tibble(value = nr / (nr + ns), parameter = "f_r")

    f_plot <- d %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x=.data$value,fill=.data$parameter)) +
      #geom_density(alpha = .3) +
      ggplot2::geom_histogram(alpha = .7, binwidth = 0.01) +
      #scale_color_manual(values = c("rho_s"=sensitive_color, "rho_r"=resistant_color)) +
      ggplot2::scale_fill_manual(values = c("f_r"=resistant_color)) +
      my_ggplot_theme() +
      ggplot2::labs(fill = "", col="", x="Resistant population fraction") +
      ggplot2::lims(x = c(-0.05,1.05))

    plots$f_plot = f_plot
  }

  w = c(1, rep(1 / (length(plots)), length(plots) - 1))

  patchwork::wrap_plots(plots, ncol = 1, heights = w)
}

# Utils

get_data_for_two_pop_plot <- function(x, alpha) {
  fit <- x$two_pop_fit

  factor_size <- x$metadata$factor_size # factor size

  # Produce ro quantiles
  rho_samples <- get_parameters(fit, par_list = c("rho_s", "rho_r"))

  rho_quantiles <- rho_samples %>%
    dplyr::group_by(.data$parameter) %>%
    dplyr::summarise(low = stats::quantile(.data$value, alpha / 2, na.rm=TRUE), mid = stats::quantile(.data$value, .5, na.rm=TRUE), high = stats::quantile(.data$value, 1 - alpha / 2, na.rm=TRUE))
  rho_s_quantiles <- rho_quantiles %>% dplyr::filter(.data$parameter == "rho_s")
  rho_r_quantiles <- rho_quantiles %>% dplyr::filter(.data$parameter == "rho_r")

  median_t0_r <- get_parameter(fit, "t0_r") %>%
    dplyr::pull(.data$value) %>%
    stats::median()

  median_t_end <- get_parameter(fit, "t_end") %>%
    dplyr::pull(.data$value) %>%
    stats::median()

  min_t <- min(median_t0_r, min(x$counts$time))
  max_t <- max(median_t_end, max(x$counts$time))

  xs <- seq(min_t, max_t, length = 2000)

  median_ns <- get_parameter(fit, "ns[1]") %>%
    dplyr::pull(.data$value) %>%
    stats::median()

  func <- two_pops_evo

  ylow <- lapply(xs, func, ns = median_ns, t0_r = median_t0_r, rho_s = -rho_s_quantiles$low, rho_r = rho_r_quantiles$low)
  ymid <- lapply(xs, func, ns = median_ns, t0_r = median_t0_r, rho_s = -rho_s_quantiles$mid, rho_r = rho_r_quantiles$mid)
  yhigh <- lapply(xs, func, ns = median_ns, t0_r = median_t0_r, rho_s = -rho_s_quantiles$high, rho_r = rho_r_quantiles$high)

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
