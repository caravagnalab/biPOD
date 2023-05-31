
my_ggplot_theme = function(cex_opt = 1)
{
  ggplot2::theme_bw(base_size = 10 * cex_opt) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(.3 * cex_opt, "cm"),
      panel.background = ggplot2::element_rect(fill = 'white')
    )
}

get_group_colors = function()
{
  color = c(
    "steelblue",
    'chocolate',
    'darkgreen',
    'firebrick3',
    'turquoise4',
    "goldenrod3"
  )
  color
}

add_shadow_to_plot = function(x, base_plot) {

  if (is.null(x$metadata$breakpoints)) {
    highlights <- dplyr::tibble(
      group = 0,
      from = min(x$counts$time),
      to = max(x$counts$time)
    )
  } else {
    ngroup <- length(x$metadata$breakpoints) + 1
    highlights <- dplyr::tibble(
      group = seq(0,ngroup-1, by=1),
      from = c(min(x$counts$time), x$metadata$breakpoints),
      to = c(x$metadata$breakpoints, max(x$counts$time))
    )
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
    ggplot2::scale_fill_manual(values = get_group_colors()) +
    ggplot2::guides(fill = ggplot2::guide_legend('', override.aes = list(alpha = 1)))

  base_plot
}

get_normalized_density <- function(values, max_value = 1) {
  # compute density and normalize it
  dens <- stats::density(values)
  dens$y <- dens$y * max_value / max(dens$y)

  # create tibble
  df <- dplyr::tibble(x=dens$x, y=dens$y)
  df
}

get_data_for_plot = function(x, alpha) {
  fit <- x$fit

  # Get fit info
  if (is.null(x$metadata$breakpoints)) {
    G <- 1
    breakpoints = array(0, dim=c(0))
  } else {
    G = length(x$metadata$breakpoints) + 1
    breakpoints <- x$metadata$breakpoints
  }

  factor_size <- x$metadata$factor_size # factor size
  t0_lower_bound <- x$metadata$t0_lower_bound

  # Plot model with t0
  if (t0_lower_bound == x$counts$time[1]) {
    median_t0 <- as.array(t0_lower_bound)
    n0 <- x$counts$count[1] / factor_size
  } else {
    median_t0 <- biPOD:::extract_parameter(x$fit, "t0") %>%
      dplyr::pull(.data$value) %>%
      stats::median()
    n0 <- 1
  }

  # Produce ro quantiles
  rho_samples <- biPOD:::extract_parameters(x$fit, par_list = paste0("rho[", 1:(length(x$metadata$breakpoints)+1), "]"))

  rho_quantiles <- rho_samples %>%
    dplyr::group_by(.data$parameter) %>%
    dplyr::summarise(low = stats::quantile(value, alpha/2), mid = stats::quantile(.data$value, .5), high = stats::quantile(value, 1 - alpha/2))

  xs <- seq(median_t0, max(x$counts$time), length=1000)

  if (x$metadata$growth_type == "logistic") {
    K <- biPOD:::extract_parameter(x$fit, "K") %>% dplyr::pull(.data$value) %>% stats::median()
    func <- biPOD:::log_growth_multiple
    ylow <- lapply(xs, func, t0=median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$low, K=K, n0=n0) %>% unlist()
    ymid <- lapply(xs, func, t0=median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$mid, K=K, n0=n0) %>% unlist()
    yhigh <- lapply(xs, func, t0=median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$high, K=K, n0=n0) %>% unlist()
  } else {
    func <- biPOD:::exp_growth
    ylow <- lapply(xs, func, t0=median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$low, n0=n0) %>% unlist()
    ymid <- lapply(xs, func, t0=median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$mid, n0=n0) %>% unlist()
    yhigh <- lapply(xs, func, t0=median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$high, n0=n0) %>% unlist()
  }

  fitted_data <- dplyr::tibble(
    x = xs,
    y = factor_size * ymid,
    ylow = factor_size * ylow,
    yhigh = factor_size * yhigh
  )

  return(fitted_data)
}

add_t0_posterior = function(base_plot, x) {

  t0_lower_bound <- x$metadata$t0_lower_bound

  # Plot model with t0
  if (t0_lower_bound == x$counts$time[1]) {

    base_plot <- base_plot +
      ggplot2::geom_segment(data = data.frame(x = t0_lower_bound, ymin = 0, ymax = max(x$counts$count)), mapping = ggplot2::aes(x=.data$x, xend=.data$x, y=.data$ymin, yend=.data$ymax), col = "darkorange")

  } else {
    values <- biPOD:::extract_parameter(x$fit, "t0") %>% dplyr::pull(.data$value)

    # Add t0 posterior
    df <- biPOD:::get_normalized_density(values, max_value = max(x$counts$count))

    base_plot <- base_plot +
      ggplot2::geom_line(data = df, mapping = ggplot2::aes(x=.data$x, y=.data$y), col = "mediumpurple") +
      ggplot2::geom_ribbon(data = df, mapping = ggplot2::aes(x=.data$x, ymin=0, ymax=.data$y), fill="mediumpurple", alpha=.5)
  }
  base_plot
}
