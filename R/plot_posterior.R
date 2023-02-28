

#' Plot growth rates posteriors
#'
#' @param x a bipod object with a 'fit' field
#' @param add_prior Boolean, indicate whether to plot also the prior distribution
#' @param labels Vector of labels for the growth rate of each group. If NULL, the standard is 'rho' plus the group index.
#'
#' @return A ggplot object containing the posterior density plots of the growth rates and the prior density plot
#' @export
#'
get_growth_rate_posteriors = function(x, add_prior = F, labels = NULL) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  # POSTERIOR
  # Obtain list of parameters to plot
  n_groups <- length(unique(x$counts$group))
  par_list <- paste0("rho[", c(1:n_groups), "]")

  d_long <- rstan::extract(x$fit, pars = par_list) %>%
    as.data.frame()

  if (is.null(labels)) {
    d_long <- d_long %>%
      `colnames<-`(par_list) %>%
      reshape2::melt(id.vars=NULL)

    d_long$variable <- factor(d_long$variable, labels = unique(c(bquote(.(d_long$variable)))))
  } else {

    if (!(length(unique(labels)) == n_groups)) stop("The number of unique labels should be equal to the number of groups")

    d_long <- d_long %>%
      `colnames<-`(labels) %>%
      reshape2::melt(id.vars=NULL)

    d_long$variable <- factor(d_long$variable, labels = unique(labels))
  }

  # plot posterior density
  # xlims <- c(min(d_long$value) * .99, max(d_long$value) * 1.01)

  # Filter data
  #d_long <- d_long %>%
  #  dplyr::filter(.data$value >= xlims[1] & .data$value <= xlims[2])

  # plot posterior density and change colors of facet box
  unique_vars <- unique(d_long$variable) %>% as.character()
  colors <- get_group_colors()

  plots <- lapply(c(1:length(unique_vars)), function(i) {
    v = unique_vars[i]

    p <- ggplot2::ggplot(d_long %>% dplyr::filter(.data$variable == v), ggplot2::aes(x=.data$value)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), alpha = .3, bins = 100) +
      ggplot2::geom_density(col = "black", linewidth = .8) +
      ggplot2::facet_wrap( ~ .data$variable, labeller = ggplot2::label_parsed) +
      my_ggplot_theme()

    p <- p + ggplot2::theme(strip.background = ggplot2::element_rect(fill = ggplot2::alpha(colors[[i]], .8)))

    return(p)
  })

  # p <- ggplot2::ggplot(d_long, ggplot2::aes(x=.data$value)) +
  #   ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), alpha = .3) +
  #   ggplot2::geom_density(col = "black", linewidth = .8) +
  #   ggplot2::facet_wrap( ~ .data$variable, labeller = ggplot2::label_parsed, scales = facet_scale)

  if (add_prior) {
    xs <- seq(xlims[1], xlims[2], length=500)
    ys <- stats::dnorm(xs)
    prior_data = data.frame(x=xs, y=ys)

    plots <- lapply(c(1:length(plots)), function(i) {
      p <- plots[[i]]

      p <- p + ggplot2::geom_line(
        data = prior_data,
        ggplot2::aes(x=.data$x, y=.data$y),
        col = "indianred3",
        size = .8)

      return(p)
    })

  }

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
    as.data.frame() %>%
    `colnames<-`(par_list) %>%
    reshape2::melt(id.vars=NULL)

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
    prior_data = data.frame(x=xs, y=ys)
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
    as.data.frame() %>%
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
    prior_data = data.frame(x=xs, y=ys)
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
