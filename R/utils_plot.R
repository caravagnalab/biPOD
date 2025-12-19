
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
    "goldenrod3",
    "mediumpurple"
  )
  color
}
