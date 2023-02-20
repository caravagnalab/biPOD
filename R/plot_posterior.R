

# plot_birth_and_death_rates_posteriors = function(x, add_prior = T) {
#   # Check input
#   if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
#   if (!("fits" %in% names(x))) stop("Input must contain a 'fits' field")
#
#   # # Check point_est and select corresponding function
#   # point_est <- match.arg(point_est)
#   # if (point_est == "median") {
#   #   m_func = stats::median
#   # } else if (point_est == "mean") {
#   #   m_func = base::mean
#   # } else {
#   #   m_func = mode.fun
#   # }
#
#   # Obtain list of parameters to plot
#   par_list = c("lambda", "mu")
#
#   # Plot for every fit
#   plots <- lapply(names(x$fits), function(fit_name) {
#     fit <- x$fits[[fit_name]]
#
#     # Prepare data
#     d_long <- rstan::extract(fit, pars = par_list) %>%
#       as.data.frame() %>%
#       dplyr::rename_at(par_list, ~paste0(par_list, gsub("fit", "", fit_name))) %>%
#       reshape2::melt(id.vars=NULL)
#
#     idx <- as.numeric(gsub("fit", "", fit_name))
#     d_long$variable <- factor(d_long$variable, ordered = TRUE, labels = c(bquote(lambda[.(idx)]), bquote(mu[.(idx)])))
#
#     if (x$fit_info$prior == "invgamma") {
#       xlims <- c(0 - 0.1, max(d_long$value) + 0.1)
#     } else {
#       xlims <- c(x$fit_info$a - 0.1, x$fit_info$b + 0.1)
#     }
#
#     # Filter data
#     d_long <- d_long %>%
#       dplyr::filter(.data$value >= xlims[1] & .data$value <= xlims[2])
#
#     # plot posterior density
#     p <- ggplot2::ggplot(d_long, ggplot2::aes(x=.data$value)) +
#       ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = .01, alpha = .3) +
#       ggplot2::geom_density(col = "forestgreen", size = .8) +
#       ggplot2::facet_wrap( ~ .data$variable, labeller = ggplot2::label_parsed)
#
#     if (add_prior) {
#       # Prepare prior data
#       prior_data <- prepare_prior_data(x, par = "lambda")
#
#       # Filter prior data
#       prior_data <- prior_data %>%
#         dplyr::filter(.data$x >= xlims[1] & .data$x <= xlims[2])
#
#       # Plot prior
#       p <- p +
#         ggplot2::geom_line(
#           data = prior_data,
#           ggplot2::aes(x=.data$x, y=.data$y),
#           col = "indianred3",
#           linewidth = .8
#         )
#     }
#
#     # Add style
#     p <- p +
#       ggplot2::labs(
#         y = 'density',
#         x = "value"
#       ) +
#       ggplot2::coord_cartesian(xlim=xlims) +
#       my_ggplot_theme()
#   })
#
#   plots <- ggpubr::ggarrange(plotlist = plots, ncol = 1)
#   plots
# }

#' Plot growth rates posteriors
#'
#' @param x a bipod object with a 'fit' field
#' @param add_prior Boolean, indicate whether to plot also the prior distribution
#'
#' @return A ggplot object containing the posterior density plots of the growth rates and the prior density plot
#' @export
#'
plot_growth_rate_posteriors = function(x, add_prior = F) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  # POSTERIOR
  # Obtain list of parameters to plot
  n_groups <- length(unique(x$counts$group))
  par_list <- paste0("rho[", c(1:n_groups), "]")

  d_long <- rstan::extract(x$fit, pars = par_list) %>%
    as.data.frame() %>%
    `colnames<-`(par_list) %>%
    reshape2::melt(id.vars=NULL)

  d_long$variable <- factor(d_long$variable, labels = unique(c(bquote(.(d_long$variable)))))

  # plot posterior density
  xlims <- c(min(d_long$value) * .99, max(d_long$value) * 1.01)
  bw = (xlims[2] - xlims[1]) / 100

  # Filter data
  d_long <- d_long %>%
    dplyr::filter(.data$value >= xlims[1] & .data$value <= xlims[2])

  # plot posterior density
  p <- ggplot2::ggplot(d_long, ggplot2::aes(x=.data$value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = bw, alpha = .3) +
    ggplot2::geom_density(col = "forestgreen", linewidth = .8) +
    ggplot2::facet_wrap( ~ .data$variable, labeller = ggplot2::label_parsed)

  if (add_prior) {
    xs <- seq(xlims[1], xlims[2], length=500)
    ys <- stats::dnorm(xs)
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
    ggplot2::coord_cartesian(xlim=xlims) +
    my_ggplot_theme()

  p

  # # Obtain list of parameters to plot
  # n_groups <- length(unique(x$counts$group))
  # par_list <- paste0("rho[", c(1:n_groups), "]")
  #
  # # Prepare data
  # d_long <- rstan::extract(x$fit, pars = par_list) %>%
  #   as.data.frame() %>%
  #   `colnames<-`(par_list) %>%
  #   reshape2::melt(id.vars=NULL)
  #
  # print(d_long)
  #
  # ####
  #
  # plots <- lapply(par_list, function(par) {
  #
  #   # Prepare data
  #   d_long <- rstan::extract(x$fit, pars = par_list) %>%
  #     as.data.frame() %>%
  #     dplyr::rename_at(par_list, ~paste0(par_list, gsub("fit", "", fit_name))) %>%
  #     reshape2::melt(id.vars=NULL)
  #
  #   idx <- as.numeric(gsub("fit", "", fit_name))
  #   d_long$variable <- factor(d_long$variable, labels = c(bquote(rho[.(idx)])))
  #
  #   # plot posterior density
  #   xlims <- c(min(d_long$value) * .99, max(d_long$value) * 1.01)
  #   bw = (xlims[2] - xlims[1]) / 100
  #
  #   # Filter data
  #   d_long <- d_long %>%
  #     dplyr::filter(.data$value >= xlims[1] & .data$value <= xlims[2])
  #
  #   # plot posterior density
  #   p <- ggplot2::ggplot(d_long, ggplot2::aes(x=.data$value)) +
  #     ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = bw, alpha = .3) +
  #     ggplot2::geom_density(col = "forestgreen", size = .8) +
  #     ggplot2::facet_wrap( ~ .data$variable, labeller = ggplot2::label_parsed)
  #
  #   if (add_prior) {
  #     # Prepare prior data
  #     prior_data <- prepare_prior_data(x, par = "ro")
  #
  #     # Filter prior data
  #     prior_data <- prior_data %>%
  #       dplyr::filter(.data$x >= xlims[1] & .data$x <= xlims[2])
  #
  #     # Create prior of ro sampling from the two priors over lambda and mu
  #     if (x$fit_info$prior == "uniform") {
  #       p <- p +
  #         ggplot2::geom_line(
  #           data = prior_data,
  #           ggplot2::aes(x=.data$x, y=.data$y),
  #           col = "indianred3",
  #           size = .8)
  #     } else if (x$fit_info$prior == "invgamma") {
  #       p <- p +
  #         ggplot2::geom_density(
  #           data = prior_data,
  #           ggplot2::aes(x=.data$x),
  #           col = "indianred3",
  #           size = .8)
  #
  #     } else {
  #       cli::cli_alert_danger("The prior {.var x$fit_info$prior} has not been recognized")
  #     }
  #   }
  #
  #   # Add style
  #   p <- p +
  #     ggplot2::labs(
  #       y = 'density',
  #       x = "value"
  #     ) +
  #     ggplot2::coord_cartesian(xlim=xlims) +
  #     my_ggplot_theme()
  # })
  #
  # if (same_scale) {
  #
  #   limits <- lapply(names(x$fits), function(fit_name) {
  #     fit <- x$fits[[fit_name]]
  #     # Prepare data
  #     d_long <- rstan::extract(fit, pars = par_list) %>%
  #       as.data.frame() %>%
  #       dplyr::rename_at(par_list, ~paste0(par_list, gsub("fit", "", fit_name))) %>%
  #       reshape2::melt(id.vars=NULL)
  #
  #     xlims <- c(min(d_long$value) * .99, max(d_long$value) * 1.01)
  #     return(xlims)
  #   })
  #
  #   min_x <- min(unlist(limits))
  #   max_x <- max(unlist(limits))
  #
  #   plots <- lapply(plots, function(p) {
  #     p <- p +
  #       ggplot2::coord_cartesian(xlim=c(min_x, max_x))
  #   })
  # }
  #
  # plots <- ggpubr::ggarrange(plotlist = plots, ncol = 1)
  # plots
}

# # Function to compute the mode
# mode.fun = function(x) {
#   print("TODO")
#   return(0)
# }
#
# # Triangular distribution
# pdf_triangular <- function(z, a, b) {
#   prob <- rep(0,length(z))
#   prob[(z >= a) & (z <= (a+b)/2)] <- 2*(z[(z >= a) & (z <= (a+b)/2)] - a)/((b-a)^2)
#   prob[(z > (a+b)/2) & (z <= b)] <- 2*(b - z[(z > (a+b)/2) & (z <= b)])/((b-a)^2)
#   return(prob)
# }
#
# prepare_prior_data = function(x, par) {
#   par_list = c(par)
#
#   if (par %in% c("ro")) {
#     limits <- lapply(x$fits, function(fit) {
#       d_long <- rstan::extract(fit, pars = par_list) %>%
#         as.data.frame() %>%
#         reshape2::melt(id.vars=NULL)
#
#       x_max <- max(d_long$value)
#       x_min <- min(d_long$value)
#       return(c(x_min, x_max))
#     })
#     x_min <- min(unlist(limits))
#     x_max <- max(unlist(limits))
#
#     if (x$fit_info$prior == "uniform") {
#       xs = seq(x_min, x_max, length=1000)
#       prior_data <- data.frame(
#         x = xs,
#         y = pdf_triangular(xs, x$fit_info$a, x$fit_info$b)
#       ) %>%
#         na.omit() %>%
#         dplyr::filter(dplyr::between(x, x_min, x_max))
#
#     } else if (x$fit_info$prior == "invgamma") {
#       lambda_prior <- invgamma::rinvgamma(10000, x$fit_info$a, x$fit_info$b)
#       mu_prior <- invgamma::rinvgamma(10000, x$fit_info$a, x$fit_info$b)
#
#       prior_data <- data.frame(x = lambda_prior - mu_prior) %>%
#         na.omit() %>%
#         dplyr::filter(dplyr::between(x, x_min, x_max))
#
#     } else {
#       cli::cli_alert_danger("The prior {.var x$fit_info$prior} has not been recognized")
#     }
#
#   } else if (par %in% c("lambda", "mu")) {
#     if (x$fit_info$prior == "uniform") {
#       xs <- seq(x$fit_info$a - 0.1, x$fit_info$b + 0.1, length=1000)
#       prior_data <- data.frame(
#         x = xs,
#         y = dunif(xs, x$fit_info$a, x$fit_info$b)
#       ) %>%
#         na.omit()
#     } else if (x$fit_info$prior == "invgamma") {
#       # plot prior, 90 % of the prior if possible, or at max twice the max of the
#       x_lim <- invgamma::qinvgamma(.9, x$fit_info$a, x$fit_info$b)
#
#       xs <- seq(0.001, x_lim, length=1000)
#       prior_data <- data.frame(
#         x = xs,
#         y = invgamma::dinvgamma(xs, x$fit_info$a, x$fit_info$b)
#       ) %>%
#         na.omit()
#     } else {
#       cli::cli_alert_danger("The prior {.var x$fit_info$prior} has not been recognized")
#     }
#   } else {
#     stop("par not recognized")
#   }
#
#   return(prior_data)
# }
#
#


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
  xlims <- c(min(d_long$value) * .99, max(d_long$value) * 1.01)
  bw = (xlims[2] - xlims[1]) / 100

  # Filter data
  d_long <- d_long %>%
    dplyr::filter(.data$value >= xlims[1] & .data$value <= xlims[2])

  # plot posterior density
  p <- ggplot2::ggplot(d_long, ggplot2::aes(x=.data$value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = bw, alpha = .3) +
    ggplot2::geom_density(col = "forestgreen", linewidth = .8) +
    ggplot2::facet_wrap( ~ .data$variable, labeller = ggplot2::label_parsed)

  if (add_prior) {
    xmin <- min(x$fit_info$t0_lower_bound, min(d_long$value)) - 3*bw
    xmax <- max(c(min(x$counts$time), max(d_long$value))) + 3*bw
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
      x = "value"
    ) +
    # ggplot2::coord_cartesian(xlim=xlims) +
    my_ggplot_theme()

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

  d_long$variable <- factor(d_long$variable, labels = c(bquote(t[0])))

  # plot posterior density
  xlims <- c(min(d_long$value) * .99, max(d_long$value) * 1.01)
  bw = (xlims[2] - xlims[1]) / 100

  # Filter data
  d_long <- d_long %>%
    dplyr::filter(.data$value >= xlims[1] & .data$value <= xlims[2])

  # plot posterior density
  p <- ggplot2::ggplot(d_long, ggplot2::aes(x=.data$value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = bw, alpha = .3) +
    ggplot2::geom_density(col = "forestgreen", linewidth = .8) +
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
