
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

  d_long <- rstan::extract(x$fit, pars = par_list) %>%
    dplyr::as_tibble()

  if (!(is.null(legend_labels))) {
    if (!(length(unique(legend_labels)) == n_groups)) stop(glue::glue("The number of unique labels should be equal to the number of groups, which is {n_groups}"))
    par_list <- unique(legend_labels)
  }

  d_long <- d_long %>%
    `colnames<-`(par_list) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "variable", values_to = "value")

  # get densities for each variable

  densities <- lapply(c(1:length(par_list)), function(i) {
    v = par_list[i]

    values <- d_long %>% dplyr::filter(.data$variable == v) %>% dplyr::select(.data$value) %>% unlist() %>% as.numeric()

    d <- get_normalized_density(values, max_value = 1) %>% dplyr::mutate(group = v)
    return(d)
  })

  densities <- do.call(rbind, densities)

  # plot each one
  colors <- get_group_colors()
  colors <- colors[1:n_groups]

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
get_growth_rate_posteriors = function(x, labels = NULL) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  # POSTERIOR
  # Obtain list of parameters to plot
  n_groups <- length(unique(x$counts$group))
  par_list <- paste0("rho[", c(1:n_groups), "]")

  if (!(is.null(labels))) {
    if (!(length(labels) == n_groups)) stop(glue::glue("The number of labels should be equal to the number of groups, which is {n_groups}"))
  } else {
    labels <- par_list
  }

  d_long <- rstan::extract(x$fit, pars = par_list) %>%
    dplyr::as_tibble() %>%
    `colnames<-`(labels) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "variable", values_to = "value")

  # plot posterior density and change colors of facet box
  unique_vars <- unique(d_long$variable) %>% as.character()
  colors <- get_group_colors()

  plots <- lapply(c(1:length(unique_vars)), function(i) {
    v = unique_vars[i]

    filtered_d <- d_long %>% dplyr::filter(.data$variable == v)

    p <- ggplot2::ggplot(filtered_d, mapping = ggplot2::aes(x = .data$value)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), alpha = .3, bins = 100) +
      ggplot2::geom_density(col = "black", linewidth = .8) +
      ggplot2::facet_wrap( ~ .data$variable, labeller = ggplot2::label_parsed) +
      my_ggplot_theme()

    p <- p + ggplot2::theme(strip.background = ggplot2::element_rect(fill = ggplot2::alpha(colors[[i]], .8)))

    return(p)
  })

  # Add style
  plots <- lapply(c(1:length(plots)), function(i) {
    p <- plots[[i]]
    p <- p +
      ggplot2::labs(
        y = 'density',
        x = ''
      )

    return(p)
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

  # Obtain list of parameters to plot
  par_list <- c("t0")

  d_long <- rstan::extract(x$fit, pars = par_list) %>%
    dplyr::as_tibble() %>%
    `colnames<-`(par_list) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "variable", values_to = "value")

  d_long$variable <- factor(d_long$variable, labels = c(bquote(t[0])))

  # plot posterior density
  #xlims <- c(min(d_long$value) * .99, max(d_long$value) * 1.01)

  # Filter data
  #d_long <- d_long %>%
  #  dplyr::filter(.data$value >= xlims[1] & .data$value <= xlims[2])

  # plot posterior density
  p <- ggplot2::ggplot(d_long, ggplot2::aes(x=.data$value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 100, alpha = .3) +
    ggplot2::geom_density(col = "black", linewidth = .8) +
    ggplot2::facet_wrap( ~ .data$variable, labeller = ggplot2::label_parsed)

  if (add_prior) {
    xmin <- min(x$fit_info$t0_lower_bound, min(d_long$value)) - 0.1
    xmax <- max(c(min(x$counts$time), max(d_long$value))) + 0.1
    xs <- seq(xmin, xmax, length=500)
    ys <- dunif(xs, x$fit_info$t0_lower_bound, min(x$counts$time))
    prior_data = dplyr::tibble(x=xs, y=ys)
    p <- p + ggplot2::geom_line(
      data = prior_data,
      ggplot2::aes(x=.data$x, y=.data$y),
      col = "indianred3",
      size = .8)
  }

  # Add style
  p <- p +
    ggplot2::labs(
      y = 'density',
      x = ''
    ) +
    # ggplot2::coord_cartesian(xlim=xlims) +
    my_ggplot_theme() +
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
  if (!(x$fit_info$growth_type == "logistic")) stop("Carrying capacity posterior is available only for 'logistic' model!")

  # Obtain list of parameters to plot
  par_list <- c("K")

  d_long <- rstan::extract(x$fit, pars = par_list) %>%
    dplyr::as_tibble() %>%
    `colnames<-`(par_list) %>%
    reshape2::melt(id.vars=NULL)

  d_long$variable <- factor(d_long$variable, labels = "K")

  # plot posterior density
  xlims <- c(min(d_long$value) * .99, max(d_long$value) * 1.01)
  bw = (xlims[2] - xlims[1]) / 100

  # Filter data
  d_long <- d_long %>%
    dplyr::filter(.data$value >= xlims[1] & .data$value <= xlims[2])

  # plot posterior density
  p <- ggplot2::ggplot(d_long, ggplot2::aes(x=.data$value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = bw, alpha = .3) +
    ggplot2::geom_density(col = "black", linewidth = .8) +
    ggplot2::facet_wrap( ~ .data$variable, labeller = ggplot2::label_parsed)

  if (add_prior) {
    prior_K <- x$fit_info$prior_K
    xmin <- prior_K * .95
    xmax <- max(d_long$value) + 3*bw
    xs <- seq(xmin, xmax, length=500)
    ys <- stats::dnorm(xs, mean=prior_K, sd = prior_K * .1)
    prior_data = dplyr::tibble(x=xs, y=ys)
    p <- p + ggplot2::geom_line(
      data = prior_data,
      ggplot2::aes(x=.data$x, y=.data$y),
      col = "indianred3",
      size = .8)
  }

  # Add style
  p <- p +
    ggplot2::labs(
      y = 'density',
      x = "value"
    ) +
    # ggplot2::coord_cartesian(xlim=xlims) +
    my_ggplot_theme()

  p
}
