# Plots of fitted growth models, from fixture draws

test_that("plot_ribbon ignores yrep draws marked -1 and draws breakpoint shadows", {
  chain <- function(shift) list(`yrep[1]` = c(100 + shift, 110 + shift, -1), `yrep[2]` = c(50, 55, 60) + shift,
                                `yrep[3]` = c(-1, -1, 20 + shift))
  x <- list(
    counts = data.frame(time = c(0, 10, 20), count = c(105, 52, 20)),
    growth_fit = list(fit = list(draws = list(chain(0), chain(1)))),
    metadata = list(breakpoints = 12)
  )
  p <- plot_ribbon(x)
  ribbon <- ggplot2::layer_data(p, 1)
  expect_true(all(ribbon$ymin >= 0))
  expect_equal(ggplot2::layer_data(p, 2)$y[3], 20.5)
})
