
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


#' Produces a list of plots, representing the posteriors, one for each growth rate.
#'
#' @param x a bipod object with a 'fit' field
#' @param labels Vector of labels for the growth rate of each group. If NULL, the standard is 'rho' plus the group index.
#'
#' @return A ggplot object containing the posterior density plots of the growth rates and the prior density plot
#' @export
#'
get_growth_rate_posteriors = function(x) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  # POSTERIOR
  # Obtain list of parameters to plot
  n_groups <- length(unique(x$counts$group))
  par_list <- paste0("rho[", c(1:n_groups), "]")
  colors <- biPOD:::get_group_colors()

  plots <- lapply(1:length(par_list), function(i) {
    par_name = par_list[i]
    color <- colors[i]

    p <- biPOD:::plot_posterior(x=x, x_fit=x$fit, par_name=par_name, color=color)
    p
  })

  return(plots)
}

#' Plot posteriors of t0
#'
#' @param x a bipod object with a 'fit' field
#' @param add_prior Boolean, indicate whether to plot also the prior distribution
#'
#' @return A ggplot object containing the posterior density plots of the growth rates and the prior density plot
#' @export
#'
plot_t0_posterior = function(x, add_prior = F) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  t0_lower_bound <- x$metadata$t0_lower_bound

  # Plot model with t0
  if (t0_lower_bound == x$counts$time[1]) {
    stop("t0 has been set manually as the first time step. Posterior distribution not available.")
  } else {
    p <- biPOD:::plot_posterior(x, x$fit, "t0", "darkorange")
  }

  if (add_prior) {
    xmin <- x$metadata$t0_lower_bound - 0.1
    xmax <- min(x$counts$time) + 0.1
    xs <- seq(xmin, xmax, length=500)
    ys <- stats::dunif(xs, x$metadata$t0_lower_bound, min(x$counts$time))
    prior_data = dplyr::tibble(x=xs, y=ys)
    p <- p + ggplot2::geom_line(
      data = prior_data,
      ggplot2::aes(x=.data$x, y=.data$y),
      col = "indianred3",
      size = .8)
  }

  # Add style
  p <- p +
    ggplot2::labs( y = 'density', x = '') +
    biPOD:::my_ggplot_theme() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = ggplot2::alpha("darkorange", .8)))

  p
}

#' Plot posteriors of carrying capacity K
#'
#' @param x a bipod object with a 'fit' field
#' @param add_prior Boolean, indicate whether to plot also the prior distribution
#'
#' @return A ggplot object containing the posterior density plots of the growth rates and the prior density plot
#' @export
#'
plot_carrying_capacity_posterior = function(x, add_prior = F) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")
  if (!(x$metadata$growth_type == "logistic")) stop("Carrying capacity posterior is available only for 'logistic' model!")

  # plot posterior density
  p <- biPOD:::plot_posterior(x, x$fit, par_name = "K", color = "seagreen")

  if (add_prior) {
    prior_K <- x$metadata$prior_K
    xmin <- prior_K * .95
    xmax <- prior_K * 10
    xs <- seq(xmin, xmax, length=500)
    ys <- stats::dnorm(xs, mean=prior_K, sd = prior_K)
    prior_data = dplyr::tibble(x=xs, y=ys)
    p <- p + ggplot2::geom_line(
      data = prior_data,
      ggplot2::aes(x=.data$x, y=.data$y),
      col = "indianred3",
      size = .8)
  }

  # Add style
  p <- p +
    ggplot2::labs( y = 'density', x = "value" ) +
    # ggplot2::coord_cartesian(xlim=xlims) +
    biPOD:::my_ggplot_theme()

  p
}
