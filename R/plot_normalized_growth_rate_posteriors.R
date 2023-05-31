
#' Produces a plot of the posteriors for the growth rates.
#' The posteriors are renormalized such that their highest peaks is at 1.
#'
#' @param x a bipod object with a 'fit' field
#' @param add_prior Boolean, indicate whether to plot also the prior distribution
#' @param legend_labels Vector containing a name for each unique fitted parameters. Default is 'rho'
#' @param legend_title Title for the legend. Default is "group"
#'
#' @return A ggplot object containing the posterior density plots of the growth rates and the prior density plot
#' @export
#'
plot_normalized_growth_rate_posteriors = function(x, add_prior = F, legend_labels = NULL, legend_title = "group") {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  # POSTERIOR
  # Obtain list of parameters to plot
  n_groups <- length(unique(x$counts$group))
  par_list <- paste0("rho[", c(1:n_groups), "]")

  d_long <- biPOD:::extract_parameters(x$fit, par_list = par_list)

  # get densities for each variable
  densities <- lapply(c(1:length(par_list)), function(i) {
    v = par_list[i]
    values <- d_long %>% dplyr::filter(.data$parameter == v) %>% dplyr::pull(.data$value)
    d <- biPOD:::get_normalized_density(values, max_value = 1) %>% dplyr::mutate(group = v)
    return(d)
  })
  densities <- do.call(rbind, densities)

  # plot each one
  colors <- biPOD:::get_group_colors()
  colors <- colors[1:n_groups]

  if (!(is.null(legend_labels))) {
    if (!(length(unique(legend_labels)) == n_groups)) stop(glue::glue("The number of unique labels should be equal to the number of groups, which is {n_groups}"))
    par_list <- unique(legend_labels)
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = densities, mapping = ggplot2::aes(x=.data$x, y=.data$y, color = .data$group)) +
    ggplot2::geom_ribbon(data = densities, mapping = ggplot2::aes(x=.data$x, y=.data$y, fill = .data$group, ymin = 0, ymax=.data$y), alpha = .3) +
    ggplot2::scale_fill_manual(values = colors, labels = par_list) +
    ggplot2::scale_color_manual(values = colors, labels = par_list) +
    ggplot2::guides(fill=ggplot2::guide_legend(title=legend_title), color=ggplot2::guide_legend(title=legend_title))

  # Add prior
  if (add_prior) {
    # plot at least between -1 and 1
    xmin <- if(min(densities$x) <= -1) min(densities$x) else -1
    xmax <- if(max(densities$x) >= 1) max(densities$x) else 1

    x <- seq(xmin, xmax, length = 500)
    y <- stats::dnorm(x)
    y_norm <- y / max(y)

    prior_data <- dplyr::tibble(x=x, y=y_norm)
    p <- p +
      ggplot2::geom_line(data = prior_data, mapping = ggplot2::aes(x=.data$x, y=.data$y), col = "darkred")
  }

  # Add style
  p <- p + my_ggplot_theme()

  return(p)
}
