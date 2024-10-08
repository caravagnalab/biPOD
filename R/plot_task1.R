#' Plot the Fit Over the Input Data
#'
#' Creates a plot showing the fit of the model over the input data.
#' The plot can include a zoomed-in view, a secondary x-axis, and highlights for different time windows.
#' Optionally, it can also display the posterior distribution for the parameter `t0`.
#'
#' @param x A `bipod` object that contains a 'fit' field.
#' @param CI Numeric value between 0 and 1 representing the confidence interval for the fit. (default is 0.95)
#' @param legend_labels Vector of labels for each unique group in `x$counts$group`. (default is NULL)
#' @param legend_title Title for the legend. (default is "group")
#' @param zoom Logical value indicating whether to present a zoomed-in view of the observations. (default is TRUE)
#' @param full_process Logical value indicating whether to include the posterior for `t0` in the plot. (default is TRUE)
#' @param sec_axis_breaks Vector of values specifying where breaks should appear on the secondary x-axis. (default is NULL)
#' @param sec_axis_labels Vector of labels for the secondary x-axis. (default is NULL)
#' @param t0_posterior_color Color to use for the posterior distribution of `t0`. (default is 'darkorange')
#' @param shadows_colors Colors to use for highlighting different time windows. (default is NULL)
#'
#' @return A `ggplot2` object displaying the fit over the input data.
#' @export
plot_fit <- function(x,
                     CI = .95,
                     legend_labels = NULL,
                     legend_title = "group",
                     zoom = TRUE,
                     full_process = TRUE,
                     sec_axis_breaks = NULL,
                     sec_axis_labels = NULL,
                     t0_posterior_color = "darkorange",
                     shadows_colors = NULL) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  alpha <- 1 - CI
  fitted_data <- get_data_for_plot(x, alpha = alpha)

  # Plot data
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(x$counts, mapping = ggplot2::aes(x = .data$time, y = .data$count)) + # original points
    ggplot2::geom_line(fitted_data, mapping = ggplot2::aes(x = .data$x, y = .data$y), col = "black") +
    ggplot2::geom_ribbon(fitted_data, mapping = ggplot2::aes(x = .data$x, y = .data$y, ymin = .data$ylow, ymax = .data$yhigh), fill = "black", alpha = .3) +
    my_ggplot_theme()

  # add highlights
  p <- add_shadow_to_plot(x, base_plot = p, colors = shadows_colors)

  # add t0 posterior
  if (full_process) p <- add_t0_posterior(base_plot = p, x = x, color = t0_posterior_color)

  # change legend
  if (!(is.null(legend_labels))) {
    if (!(length(legend_labels) == length(unique(x$counts$group)))) stop("The number of labels for the legend must be equal to the number of groups")
    p <- p + ggplot2::scale_fill_manual(values = get_group_colors(), labels = legend_labels)
  }
  p <- p + ggplot2::guides(fill = ggplot2::guide_legend(title = legend_title, override.aes = list(alpha = 1)))

  # Add zoom
  if (zoom) {
    zoom_limits <- c(min(x$counts$time), max(x$counts$time))
    p <- p + ggforce::facet_zoom(xlim = zoom_limits)
  }

  # Add secondary axis
  if (xor(is.null(sec_axis_breaks), is.null(sec_axis_labels))) stop("To plot secondary x axis both sec_axis_breaks and sex_axis_labels are needed")
  if (!(is.null(sec_axis_breaks))) {
    p <- p +
      ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., breaks = sec_axis_breaks, labels = sec_axis_labels)) +
      ggplot2::theme(axis.text.x.top = ggplot2::element_text(angle = 90, vjust = 0.5))
  }

  return(p)
}

#' Plot Normalized Posterior Distributions of Growth Rates
#'
#' Creates a plot of the normalized posterior distributions of the growth rates.
#' The distributions are renormalized so that their highest peak is at 1. Optionally, the prior distribution can be included in the plot.
#'
#' @param x A `bipod` object that contains a 'fit' field.
#' @param add_prior Logical value indicating whether to plot the prior distribution along with the posterior distributions. (default is FALSE)
#' @param legend_labels Vector of labels for each unique fitted parameter. (default is NULL, in which case 'rho' is used)
#' @param legend_title Title for the legend. (default is 'group')
#' @param colors Colors to use for the different growth rate posteriors. (default is NULL, which uses the default group colors)
#'
#' @return A `ggplot2` object containing the posterior density plots of the growth rates and optionally the prior density plot.
#' @export
plot_normalized_growth_rate_posteriors <- function(x,
                                                   add_prior = F,
                                                   legend_labels = NULL,
                                                   legend_title = "group",
                                                   colors = NULL) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  # POSTERIOR
  # Obtain list of parameters to plot
  n_groups <- length(unique(x$counts$group))
  par_list <- paste0("rho[", c(1:n_groups), "]")

  d_long <- get_parameters(x$fit, par_list = par_list, x$metadata$sampling == "variational")

  # get densities for each variable
  densities <- lapply(c(1:length(par_list)), function(i) {
    v <- par_list[i]
    values <- d_long %>%
      dplyr::filter(.data$parameter == v) %>%
      dplyr::pull(.data$value)
    d <- get_normalized_density(values, max_value = 1) %>% dplyr::mutate(group = v)
    return(d)
  })
  densities <- do.call(rbind, densities)

  # plot each one
  if (is.null(colors)) {
    colors <- get_group_colors()
    colors <- colors[1:n_groups]
  }

  if (!(is.null(legend_labels))) {
    if (!(length(unique(legend_labels)) == n_groups)) stop(glue::glue("The number of unique labels should be equal to the number of groups, which is {n_groups}"))
    par_list <- unique(legend_labels)
  }

  x_lim <- max(abs(min(densities$x)), abs(max(densities$x)))

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = densities, mapping = ggplot2::aes(x = .data$x, y = .data$y, color = .data$group)) +
    ggplot2::geom_ribbon(data = densities, mapping = ggplot2::aes(x = .data$x, y = .data$y, fill = .data$group, ymin = 0, ymax = .data$y), alpha = .3) +
    ggplot2::scale_fill_manual(values = colors, labels = par_list) +
    ggplot2::scale_color_manual(values = colors, labels = par_list) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = legend_title), color = ggplot2::guide_legend(title = legend_title)) +
    ggplot2::xlim(-x_lim, x_lim)

  # Add prior
  if (add_prior) {
    # plot at least between -1 and 1
    xmin <- if (min(densities$x) <= -1) min(densities$x) else -1
    xmax <- if (max(densities$x) >= 1) max(densities$x) else 1

    x <- seq(xmin, xmax, length = 500)
    y <- stats::dnorm(x)
    y_norm <- y / max(y)

    prior_data <- dplyr::tibble(x = x, y = y_norm)
    p <- p +
      ggplot2::geom_line(data = prior_data, mapping = ggplot2::aes(x = .data$x, y = .data$y), col = "darkred")
  }

  # Add style
  p <- p + my_ggplot_theme()

  return(p)
}

#' Plot Posteriors of t0
#'
#' Creates a plot of the posterior distribution of the parameter `t0` from a `bipod` object. Optionally, the prior distribution can be included in the plot.
#'
#' @param x A `bipod` object that contains a 'fit' field.
#' @param add_prior Logical value indicating whether to plot the prior distribution along with the posterior distribution. (default is FALSE)
#'
#' @return A `ggplot2` object containing the posterior density plot of `t0`, and optionally the prior density plot.
#' @export
plot_t0_posterior <- function(x, add_prior = F) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  # Plot model with t0
  if (x$metadata$t0_inferred) {
    p <- plot_posterior(x, x$fit, "t0", "darkorange")
  } else {
    stop("t0 has not been inferred. Posterior distribution not available.")
  }

  if (add_prior) {
    xmin <- 100
    xmax <- min(x$counts$time) + 0.1
    xs <- seq(xmin, xmax, length = 500)
    ys <- stats::dnorm(xs, min(x$counts$time), 100)
    prior_data <- dplyr::tibble(x = xs, y = ys)
    p <- p + ggplot2::geom_line(
      data = prior_data,
      ggplot2::aes(x = .data$x, y = .data$y),
      col = "indianred3",
      size = .8
    )
  }

  # Add style
  p <- p +
    ggplot2::labs(y = "density", x = "") +
    my_ggplot_theme() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = ggplot2::alpha("darkorange", .8)))

  p
}

# Utils

get_data_for_plot <- function(x, alpha) {
  fit <- x$fit

  # Get fit info
  if (is.null(x$metadata$breakpoints)) {
    G <- 1
    breakpoints <- array(0, dim = c(0))
  } else {
    G <- length(x$metadata$breakpoints) + 1
    breakpoints <- x$metadata$breakpoints
  }

  factor_size <- x$metadata$factor_size # factor size

  # Plot model with t0
  if (!(x$metadata$t0_inferred)) {
    median_t0 <- min(x$counts$time)
    n0 <- get_parameter(x$fit, "n0", variational = x$metadata$sampling == "variational") %>%
      dplyr::pull(.data$value) %>%
      stats::median()
    n0 <- n0
  } else {
    median_t0 <- get_parameter(x$fit, "t0", x$metadata$sampling == "variational") %>%
      dplyr::pull(.data$value) %>%
      stats::median()
    n0 <- 1
  }

  # Produce ro quantiles
  rho_samples <- get_parameters(x$fit, par_list = paste0("rho[", 1:(length(x$metadata$breakpoints) + 1), "]"), x$metadata$sampling == "variational")

  rho_quantiles <- rho_samples %>%
    dplyr::group_by(.data$parameter) %>%
    dplyr::summarise(low = stats::quantile(.data$value, alpha / 2), mid = stats::quantile(.data$value, .5), high = stats::quantile(.data$value, 1 - alpha / 2))

  xs <- seq(median_t0, max(x$counts$time), length = 1000)

  if (x$metadata$growth_type == "logistic") {
    K <- get_parameter(x$fit, "K", x$metadata$sampling == "variational") %>%
      dplyr::pull(.data$value) %>%
      stats::median()
    func <- log_growth_multiple
    ylow <- lapply(xs, func, t0 = median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$low, K = K, n0 = n0) %>% unlist()
    ymid <- lapply(xs, func, t0 = median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$mid, K = K, n0 = n0) %>% unlist()
    yhigh <- lapply(xs, func, t0 = median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$high, K = K, n0 = n0) %>% unlist()
  } else {
    func <- exp_growth
    ylow <- lapply(xs, func, t0 = median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$low, n0 = n0) %>% unlist()
    ymid <- lapply(xs, func, t0 = median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$mid, n0 = n0) %>% unlist()
    yhigh <- lapply(xs, func, t0 = median_t0, t_array = as.array(breakpoints), rho_array = rho_quantiles$high, n0 = n0) %>% unlist()
  }

  fitted_data <- dplyr::tibble(
    x = xs,
    y = factor_size * ymid,
    ylow = factor_size * ylow,
    yhigh = factor_size * yhigh
  )

  return(fitted_data)
}
