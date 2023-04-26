
my_ggplot_theme = function(cex_opt = 1)
{
  ggplot2::theme_light(base_size = 10 * cex_opt) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(.3 * cex_opt, "cm"),
      panel.background = ggplot2::element_rect(fill = 'white')
    )
}

get_group_colors = function()
{
  color = c(
    "forestgreen",
    "steelblue",
    'darkblue',
    'turquoise4',
    'orange',
    'firebrick3'
  )
  color
}

add_shadow_to_plot = function(x, base_plot) {



  if (is.null(x$breakpoints)) {
    highlights <- dplyr::tibble(
      group = 0,
      from = min(x$counts$time),
      to = max(x$counts$time)
    )
  } else {
    ngroup <- length(x$breakpoints) + 1
    highlights <- dplyr::tibble(
      group = seq(0,ngroup-1, by=1),
      from = c(min(x$counts$time), x$breakpoints),
      to = c(x$breakpoints, max(x$counts$time))
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
