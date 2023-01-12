
#' Plot birth and death rates posteriors
#'
#' @param x a bipod object with a 'fit' field
#' @param point_est Character string indicating which point estimate to use. Available options are "mean", "median", and "none". Default is "mean".
#'
#' @return A ggplot object containing the posterior density plots of the birth and death rates and the prior density plot
#'
#' @export
#'
plot_birth_and_death_rates_posteriors = function(x, point_est = c("mean", "median", "none")) {
  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  assertthat::assert_that("fit" %in% names(x), msg = "Input must contain a 'fit' field")

  # Check point_est and select corresponding function
  point_est <- match.arg(point_est)
  if (point_est == "median") {
    m_func = median
  } else if (point_est == "mean") {
    m_func = mean
  } else {
    m_func = mode.fun
  }

  # Obtain list of parameters to plot
  par_list = c()
  for (p in c("lambda", "mu")) {
    par_list <- c(par_list, names(x$fit)[grep(p, names(x$fit))])
  }
  assertthat::assert_that(length(par_list) > 0, msg = "This fit does not contain parameters named 'lambda' nor 'mu")

  # Prepare data
  d_long <- rstan::extract(x$fit, pars = par_list) %>%
    as.data.frame() %>%
    reshape2::melt(id.vars=NULL)

  d <- d_long %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      point_est = point_est,
      # l = quantile(.data$value, probs[1]),
      m  = m_func(.data$value),
      # h = quantile(.data$value, probs[2])
    )

  # plot posterior density
  p <- ggplot2::ggplot(d_long, aes(x=value)) +
    ggplot2::geom_histogram(aes(y = after_stat(density)), binwidth = 0.01, alpha = .3) +
    ggplot2::geom_density(col = "forestgreen", size = .8) +
    ggplot2::facet_wrap( ~ variable)

  # Plot estimates, if requested
  if (point_est != "none") {
    p <- p +
      ggplot2::geom_vline(
        data = d,
        ggplot2::aes(xintercept = m),
        col = "forestgreen",
        linetype = 'longdash',
        size = .5,
        show.legend = FALSE
      )
      # ggplot2::annotate(
      #   geom = 'rect',
      #   xmin = d$l,
      #   xmax = d$h,
      #   ymin = 0,
      #   ymax = Inf,
      #   color = NA,
      #   fill = 'forestgreen',
      #   alpha = .4
      # ) +
      # ggplot2::facet_wrap( ~ variable)
  }

  # Plot prior
  if (x$fit_info$prior == "uniform") {
    xs <- seq(x$fit_info$a - 0.1, x$fit_info$b + 0.1, length=1000)
    prior_data <- data.frame(
      x = xs,
      y = dunif(xs, x$fit_info$a, x$fit_info$b)
    ) %>%
      na.omit()
  } else if (x$fit_info$prior == "invgamma") {
    # plot prior, 90 % of the prior if possible, or at max twice the max of the
    x_lim <- max(d_long$value)
    x_lim_prior <- invgamma::qinvgamma(.9, x$fit_info$a, x$fit_info$b)
    x_lim = if(x_lim_prior > 2 * x_lim) 2 * x_lim else x_lim_prior

    xs <- seq(0.001, x_lim, length=1000)
    prior_data <- data.frame(
      x = xs,
      y = invgamma::dinvgamma(xs, x$fit_info$a, x$fit_info$b)
    ) %>%
      na.omit()
  } else {
    cli::cli_alert_danger("The prior {.var x$fit_info$prior} has not been recognized")
  }

  p <- p +
    ggplot2::geom_line(
      data = prior_data,
      ggplot2::aes(x=x, y=y),
      col = "indianred3",
      size = .8
    )

  # Add style
  p <- p +
    ggplot2::labs(
      y = 'density',
      x = "value"
    ) +
    my_ggplot_theme()

  return(p)
}

#' Plot growth rates posteriors
#'
#' @param x a bipod object with a 'fit' field
#' @param point_est Character string indicating which point estimate to use. Available options are mean", "median", and "none". Default is "mean".
#'
#' @return A ggplot object containing the posterior density plots of the growth rates and the prior density plot
#'
#' @export
#'
plot_growth_rate_posteriors = function(x, point_est = c("mean", "median", "none")) {
  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  assertthat::assert_that("fit" %in% names(x), msg = "Input must contain a 'fit' field")

  # Check point_est and select corresponding function
  point_est <- match.arg(point_est)
  if (point_est == "median") {
    m_func = median
  } else if (point_est == "mean") {
    m_func = mean
  } else {
    m_func = mode.fun
  }

  # Obtain list of parameters to plot
  par_list = c()
  for (p in c("ro")) {
    par_list <- c(par_list, names(x$fit)[grep(p, names(x$fit))])
  }
  assertthat::assert_that(length(par_list) > 0, msg = "This fit does not contain parameters named 'ro'")

  # Prepare data
  d_long <- rstan::extract(x$fit, pars = par_list) %>%
    as.data.frame() %>%
    reshape2::melt(id.vars=NULL)

  d <- d_long %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      point_est = point_est,
      # l = quantile(.data$value, probs[1]),
      m  = m_func(.data$value),
      # h = quantile(.data$value, probs[2])
    )

  # plot posterior density
  x_max <- max(d_long$value) + 0.1
  x_min <- min(d_long$value) - 0.1
  bw = (x_max - x_min) / 100

  d_long %>%
    dplyr::filter(dplyr::between(value, x_min, x_max)) %>%
    ggplot2::ggplot(aes(x=value)) +
    ggplot2::geom_histogram(aes(y = after_stat(density)), binwidth = bw, alpha = .3) +
    ggplot2::geom_density(col = "forestgreen", size = .8) +
    ggplot2::facet_wrap( ~ variable) -> p

  # Plot estimates, if requested
  if (point_est != "none") {
    p <- p +
      ggplot2::geom_vline(
        data = d,
        ggplot2::aes(xintercept = m),
        col = "forestgreen",
        linetype = 'longdash',
        size = .5,
        show.legend = FALSE
      ) +
      # ggplot2::annotate(
      #   geom = 'rect',
      #   xmin = d$l,
      #   xmax = d$h,
      #   ymin = 0,
      #   ymax = Inf,
      #   color = NA,
      #   fill = 'forestgreen',
      #   alpha = .4
      # ) +
      ggplot2::facet_wrap( ~ variable)
  }

  # Create prior of ro sampling from the two priors over lambda and mu

  if (x$fit_info$prior == "uniform") {
    xs = seq(x_min, x_max, length=1000)
    prior_data <- data.frame(
      x = xs,
      y = pdf_triangular(xs, x$fit_info$a, x$fit_info$b)
    ) %>%
      na.omit() %>%
      dplyr::filter(dplyr::between(x, x_min, x_max))

    p <- p +
      ggplot2::geom_line(
        data = prior_data,
        ggplot2::aes(x=x, y=y),
        col = "indianred3",
        size = .8) +
      xlim(x_min, x_max)

  } else if (x$fit_info$prior == "invgamma") {
    # plot prior, 90 % of the prior if possible, or at max twice the max of the
    lambda_prior <- invgamma::rinvgamma(10000, x$fit_info$a, x$fit_info$b)
    mu_prior <- invgamma::rinvgamma(10000, x$fit_info$a, x$fit_info$b)
    prior_data <- data.frame(x = lambda_prior - mu_prior) %>%
      na.omit() %>%
      dplyr::filter(dplyr::between(x, x_min, x_max))

    p <- p +
      ggplot2::geom_density(
        data = prior_data,
        ggplot2::aes(x=x),
        col = "indianred3",
        size = .8) +
      ggplot2::scale_x_continuous(limits = c(x_min, x_max), oob = scales::squish)

  } else {
    cli::cli_alert_danger("The prior {.var x$fit_info$prior} has not been recognized")
  }

  # Add style
  p <- p +
    ggplot2::labs(
      y = 'density',
      x = "value"
    ) +
    my_ggplot_theme()

  return(p)
}

# Function to compute the mode
mode.fun = function(x) {
  print("TODO")
  return(0)
}

# Triangular distribution
pdf_triangular <- function(z, a, b) {
  prob <- rep(0,length(z))
  prob[(z >= a) & (z <= (a+b)/2)] <- 2*(z[(z >= a) & (z <= (a+b)/2)] - a)/((b-a)^2)
  prob[(z > (a+b)/2) & (z <= b)] <- 2*(b - z[(z > (a+b)/2) & (z <= b)])/((b-a)^2)
  return(prob)
}

