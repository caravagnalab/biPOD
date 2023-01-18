
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
