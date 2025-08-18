
get_group_colors <- function() {
  color <- c(
    "steelblue",
    "chocolate",
    "darkgreen",
    "firebrick3",
    "turquoise4",
    "goldenrod3",
    "mediumpurple"
  )
  color
}

add_breakpoint_shadows = function(p, shadow_breakpoints, colors = NULL) {
  shadow_breakpoints = sort(shadow_breakpoints)

  xlims = ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range

  ngroup <- length(shadow_breakpoints) + 1
  highlights <- dplyr::tibble(
    group = seq(0, ngroup - 1, by = 1),
    from = c(min(xlims), shadow_breakpoints),
    to = c(shadow_breakpoints, max(xlims))
  )

  if (is.null(colors)) {
    colors <- get_group_colors()
  }

  g <- as.character(highlights$group)
  p <- p +
    ggplot2::geom_rect(
      data = highlights,
      ggplot2::aes(
        xmin = .data$from,
        xmax = .data$to,
        ymin = 0,
        ymax = Inf,
        fill = factor(.data$group, levels = g)
      ),
      alpha = .2,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::guides(fill = ggplot2::guide_legend("", override.aes = list(alpha = 1)))

  p
}
