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
