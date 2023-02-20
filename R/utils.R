
group_contiguous <- function(x) {
  rle_x <- rle(x)
  return(rep(seq_along(rle_x$lengths) - 1, times = rle_x$lengths))
}

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
    '0' = ggplot2::alpha('forestgreen', .8),
    '1' = 'steelblue',
    '2' = 'darkblue',
    '3' = 'turquoise4',
    '4' = ggplot2::alpha('orange', .8),
    '5' = 'firebrick3'
  )
  color
}

add_shadow_to_plot = function(x, base_plot) {
  times <- (x$counts %>% dplyr::group_by(.data$group) %>% dplyr::summarise(times = max(.data$time)))$times

  ngroup <- length(times)
  n_lower <- ngroup - 1

  highlights <- data.frame(
    group = seq(0,ngroup-2, by=1),
    from = times[1:n_lower],
    to = times[2:ngroup]
  )

  g <- as.character(highlights$group)
  base_plot <- base_plot +
    ggplot2::geom_rect(
      data = highlights,
      ggplot2::aes(
        xmin = from,
        xmax = to,
        ymin = -Inf,
        ymax = Inf,
        fill = factor(group, levels = g)
      ),
      alpha = .2
    ) +
    ggplot2::scale_fill_manual(values = biPOD:::get_group_colors()) +
    ggplot2::guides(fill = ggplot2::guide_legend('', override.aes = list(alpha = 1)))

  base_plot
}
