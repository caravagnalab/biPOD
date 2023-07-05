my_ggplot_theme <- function(cex_opt = 1) {
  ggplot2::theme_bw(base_size = 10 * cex_opt) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(.3 * cex_opt, "cm"),
      panel.background = ggplot2::element_rect(fill = "white")
    )
}

get_group_colors <- function() {
  color <- c(
    "steelblue",
    "chocolate",
    "darkgreen",
    "firebrick3",
    "turquoise4",
    "goldenrod3"
  )
  color
}

add_shadow_to_plot <- function(x, base_plot, colors) {
  if (is.null(x$metadata$breakpoints)) {
    highlights <- dplyr::tibble(
      group = 0,
      from = min(x$counts$time),
      to = max(x$counts$time)
    )
  } else {
    ngroup <- length(x$metadata$breakpoints) + 1
    highlights <- dplyr::tibble(
      group = seq(0, ngroup - 1, by = 1),
      from = c(min(x$counts$time), x$metadata$breakpoints),
      to = c(x$metadata$breakpoints, max(x$counts$time))
    )
  }

  if (is.null(colors)) {
    colors <- get_group_colors()
  }

  g <- as.character(highlights$group)
  base_plot <- base_plot +
    ggplot2::geom_rect(
      data = highlights,
      ggplot2::aes(
        xmin = .data$from,
        xmax = .data$to,
        ymin = 0,
        ymax = Inf,
        fill = factor(.data$group, levels = g)
      ),
      alpha = .2
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::guides(fill = ggplot2::guide_legend("", override.aes = list(alpha = 1)))

  base_plot
}

get_normalized_density <- function(values, max_value = 1) {
  # compute density and normalize it
  dens <- stats::density(values)
  dens$y <- dens$y * max_value / max(dens$y)

  # create tibble
  df <- dplyr::tibble(x = dens$x, y = dens$y)
  df
}

get_data_for_plot <- function(x, alpha) {
  fit <- x$fit

  # Get fit info
  if (is.null(x$metadata$breakpoints)) {
    G <- 1
    breakpoints <- array(0, dim = c(0))
  } else {
    G <- length(x$metadata$breakpoints) + 1
    breakpoints <- x$metadata$breakpoints
  }

  factor_size <- x$metadata$factor_size # factor size
  t0_lower_bound <- x$metadata$t0_lower_bound

  # Plot model with t0
  if (t0_lower_bound == x$counts$time[1]) {
    median_t0 <- as.array(t0_lower_bound)
    n0 <- x$counts$count[1] / factor_size
  } else {
    median_t0 <- get_parameter(x$fit, "t0") %>%
      dplyr::pull(.data$value) %>%
      stats::median()
    n0 <- 1
  }

  # Produce ro quantiles
  rho_samples <- get_parameters(x$fit, par_list = paste0("rho[", 1:(length(x$metadata$breakpoints) + 1), "]"))

  rho_quantiles <- rho_samples %>%
    dplyr::group_by(.data$parameter) %>%
    dplyr::summarise(low = stats::quantile(.data$value, alpha / 2), mid = stats::quantile(.data$value, .5), high = stats::quantile(.data$value, 1 - alpha / 2))

  xs <- seq(median_t0, max(x$counts$time), length = 1000)

  if (x$metadata$growth_type == "logistic") {
    K <- get_parameter(x$fit, "K") %>%
      dplyr::pull(.data$value) %>%
      stats::median()
    func <- log_growth_multiple
    ylow <- lapply(xs, func, t0 = median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$low, K = K, n0 = n0) %>% unlist()
    ymid <- lapply(xs, func, t0 = median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$mid, K = K, n0 = n0) %>% unlist()
    yhigh <- lapply(xs, func, t0 = median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$high, K = K, n0 = n0) %>% unlist()
  } else {
    func <- exp_growth
    ylow <- lapply(xs, func, t0 = median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$low, n0 = n0) %>% unlist()
    ymid <- lapply(xs, func, t0 = median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$mid, n0 = n0) %>% unlist()
    yhigh <- lapply(xs, func, t0 = median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$high, n0 = n0) %>% unlist()
  }

  fitted_data <- dplyr::tibble(
    x = xs,
    y = factor_size * ymid,
    ylow = factor_size * ylow,
    yhigh = factor_size * yhigh
  )

  return(fitted_data)
}

add_t0_posterior <- function(base_plot, x, color) {
  t0_lower_bound <- x$metadata$t0_lower_bound

  # Plot model with t0
  if (t0_lower_bound == x$counts$time[1]) {
    return(base_plot)
  } else {
    values <- get_parameter(x$fit, "t0") %>% dplyr::pull(.data$value)

    # Add t0 posterior
    df <- get_normalized_density(values, max_value = max(x$counts$count))

    base_plot <- base_plot +
      ggplot2::geom_line(data = df, mapping = ggplot2::aes(x = .data$x, y = .data$y), col = color) +
      ggplot2::geom_ribbon(data = df, mapping = ggplot2::aes(x = .data$x, ymin = 0, ymax = .data$y), fill = color, alpha = .5)
  }
}

add_posterior <- function(base_plot, param, x, x_fit, color) {
  values <- get_parameter(x_fit, param) %>% dplyr::pull(.data$value)

  # Add t0 posterior
  df <- get_normalized_density(values, max_value = max(x$counts$count))

  base_plot <- base_plot +
    ggplot2::geom_line(data = df, mapping = ggplot2::aes(x = .data$x, y = .data$y), col = color) +
    ggplot2::geom_ribbon(data = df, mapping = ggplot2::aes(x = .data$x, ymin = 0, ymax = .data$y), fill = color, alpha = .5)
}

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
