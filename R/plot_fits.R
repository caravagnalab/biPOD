
#' Plot the fit over the input data.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param final_time The final time for which the fit need to be plotted.
#' @param add_title Boolean, indicating whether the plot should have a title
#' @param add_highlights Boolean, indicating whether the groups should be highlighted
#' @param plot_t0 Boolean, indicating whether to plot t0 histogram and whole evolution plot
#'
#' @returns A plot of the fit over the input data.
#' @export
plot_fit = function(x, final_time = NULL, add_title = F, add_highlights = F, plot_t0 = F) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  growth_type <- x$fit_info$growth_type

  if (growth_type == "exponential") {
    p <- plot_exponential_fit(x, final_time, add_title, plot_t0)
  } else {
    p <- plot_logistic_fit(x, final_time, add_title, plot_t0)
  }

  if (add_highlights) {
    p <- biPOD:::add_shadow_to_plot(x, base_plot = p)
  }

  return(p)
}

plot_exponential_fit = function(x, final_time, add_title, plot_t0) {
  fit <- x$fit

  # Get fit info
  G <- length(unique(x$counts$group))

  if (G == 1) {
    t_array = array(0, dim=c(0))
  } else {
    n <- G - 1
    t_array <- (x$counts %>% dplyr::group_by(.data$group) %>% dplyr::slice_tail(n=1) %>% dplyr::ungroup() %>%  dplyr::select(.data$time))$time
    t_array <- t_array[1:n]
  }

  factor_size <- x$fit_info$factor_size # factor size
  if(is.null(final_time)) final_time = max(x$counts$time) # final time

  # Plot model with t0
  t0 <- rstan::extract(fit, pars=c('t0')) %>% as.list() %>% unlist() %>% as.numeric()
  t0 <- round(t0, 2)

  rho <- rstan::extract(fit, pars=c("rho")) %>% as.data.frame()
  ro_quantiles <- apply(rho, 2, function(x) stats::quantile(x, c(0.05, 0.5, 0.95))) %>% as.data.frame() %>% t() %>% as.data.frame()

  mode_t0 <- median(t0)
  # real_t0 <- (log(factor_size) - ro_quantiles[2]*mode_t0) / (-ro_quantiles[2])
  xs <- seq(mode_t0, final_time, length=1000)
  ylow <- lapply(xs, exp_growth, t0=mode_t0, t_array = as.array(t_array), rho_array = ro_quantiles$`5%`) %>% unlist()
  ymid <- lapply(xs, exp_growth, t0=mode_t0, t_array = as.array(t_array), rho_array = ro_quantiles$`50%`) %>% unlist()
  yhigh <- lapply(xs, exp_growth, t0=mode_t0, t_array = as.array(t_array), rho_array = ro_quantiles$`95%`) %>% unlist()

  fitted_data <- data.frame(
    x = xs,
    y = factor_size * ymid,
    ylow = factor_size * ylow,
    yhigh = factor_size * yhigh
  )


  # Plot them
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(x$counts, mapping = aes(x=.data$time, y=.data$count)) + #original points
    ggplot2::geom_line(fitted_data, mapping=aes(x=.data$x, y=.data$y), col="darkgreen") +
    ggplot2::geom_ribbon(fitted_data, mapping=aes(x=.data$x, y=.data$y, ymin=.data$ylow, ymax=.data$yhigh), fill="darkgreen", alpha=.3) +
    my_ggplot_theme()

  if (plot_t0) p <- p + ggplot2::geom_histogram(data.frame(x=t0), mapping=aes(x=.data$x), binwidth = .2, fill="darkorange", alpha=.5)
  if (!plot_t0) p <- p + ggplot2::scale_x_continuous(limits = c(min(x$counts$time), final_time))

  p
}

plot_logistic_fit = function(x, final_time, add_title, plot_t0) {
  fit <- x$fit
  # Get fit info
  G <- length(unique(x$counts$group))

  if (G == 1) {
    t_array = array(0, dim=c(0))
  } else {
    n <- G - 1
    t_array <- (x$counts %>% dplyr::group_by(.data$group) %>% dplyr::slice_tail(n=1) %>% dplyr::ungroup() %>% dplyr::select(.data$time))$time
    t_array <- t_array[1:n]
  }

  factor_size <- x$fit_info$factor_size # factor size
  if(is.null(final_time)) final_time = max(x$counts$time) # final time

  # Plot model with t0
  t0 <- rstan::extract(fit, pars=c('t0')) %>% as.list() %>% unlist() %>% as.numeric()
  t0 <- round(t0, 2)

  K <- mean(rstan::extract(fit, pars=c('K')) %>% as.list() %>% unlist())

  rho <- rstan::extract(fit, pars=c("rho")) %>% as.data.frame()
  ro_quantiles <- apply(rho, 2, function(x) stats::quantile(x, c(0.05, 0.5, 0.95))) %>% as.data.frame() %>% t() %>% as.data.frame()

  mode_t0 <- median(t0)
  # real_t0 <- (log(factor_size) - ro_quantiles[2]*mode_t0) / (-ro_quantiles[2])
  xs <- seq(mode_t0, max(x$counts$time), length=1000)

  ylow <- lapply(xs, log_growth_multiple, t0=mode_t0, t_array = as.array(t_array), rho_array = ro_quantiles$`5%`, K=K) %>% unlist()
  ymid <- lapply(xs, log_growth_multiple, t0=mode_t0, t_array = as.array(t_array), rho_array = ro_quantiles$`50%`, K=K) %>% unlist()
  yhigh <- lapply(xs, log_growth_multiple, t0=mode_t0, t_array = as.array(t_array), rho_array = ro_quantiles$`95%`, K) %>% unlist()

  fitted_data <- data.frame(
    x = xs,
    y = factor_size * ymid,
    ylow = factor_size * ylow,
    yhigh = factor_size * yhigh
  )

  # Plot them
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(x$counts, mapping = aes(x=.data$time, y=.data$count)) + #original points
    #geom_line(N_rep_quantiles, mapping=aes(x=time, y=mid), col="steelblue") +
    #geom_ribbon(N_rep_quantiles, mapping=aes(x=time, y=mid, ymin=low, ymax=high), fill="steelblue", alpha=.3) +
    ggplot2::geom_histogram(data.frame(x=t0), mapping=aes(x=.data$x), binwidth = .2, fill="darkorange", alpha=.5) +
    ggplot2::geom_line(fitted_data, mapping=aes(x=.data$x, y=.data$y), col="darkgreen") +
    ggplot2::geom_ribbon(fitted_data, mapping=aes(x=.data$x, y=.data$y, ymin=.data$ylow, ymax=.data$yhigh), fill="darkgreen", alpha=.3) +
    my_ggplot_theme()

  if (plot_t0) p <- p + ggplot2::geom_histogram(data.frame(x=t0), mapping=aes(x=.data$x), binwidth = .2, fill="darkorange", alpha=.5)
  if (!plot_t0) p <- p + ggplot2::scale_x_continuous(limits = c(min(x$counts$time), final_time))

  p

}

# plot_exponential_fit_old = function(x, final_time, add_title) {
#
#   p <- ggplot2::ggplot()
#   data <- data.frame()
#   groups <- unique(x$counts$group)
#
#   for (i in 2:length(groups)) {
#     previous <- x$counts %>%
#       dplyr::filter(.data$group == groups[i-1])
#
#     previous_t <- previous$time[nrow(previous)]
#     previous_n <- previous$count[nrow(previous)]
#
#     current <- x$counts %>%
#       dplyr::filter(.data$group == groups[i])
#
#     if (i == length(groups)) {
#       final_t <- if (is.null(final_time)) current$time[nrow(current)] else final_time
#     } else {
#       final_t <- current$time[nrow(current)]
#     }
#
#     xs <- seq(0, final_t - previous_t, length = 100)
#
#     # Extract fit info
#     fit_name <- paste0("fit", groups[i])
#     fit <- x$fits[[fit_name]]
#     ros <- unname(stats::quantile(rstan::extract(fit, pars="ro")$ro, c(.05, .5, .95)))
#
#     N_low = previous_n * exp(ros[1] * xs)
#     N_medium = previous_n * exp(ros[2] * xs)
#     N_high = previous_n * exp(ros[3] * xs)
#
#     d <- data.frame(x=xs + previous_t, yl=N_low, ym=N_medium, yh=N_high)
#     p <- p +
#       ggplot2::geom_line(d, mapping=ggplot2::aes(x=.data$x, y=.data$ym), col="forestgreen") +
#       ggplot2::geom_ribbon(d, mapping=ggplot2::aes(x=.data$x, y=.data$ym, ymin=.data$yl, ymax=.data$yh), fill="forestgreen", alpha=.5)
#   }
#
#   p <- p +
#     ggplot2::geom_point(x$counts, mapping=ggplot2::aes(x=.data$time, y=.data$count)) +
#     my_ggplot_theme()
#
#   if (add_title) {
#     p <- p + ggplot2::labs(
#       title = paste("Exponential fit", x$sample),
#       x = "time",
#       y = "count"
#     )
#   } else {
#     p <- p + ggplot2::labs(
#       x = "time",
#       y = "count"
#     )
#   }
#
#   p
# }

plot_logistic_fit_old = function(x, final_time, add_title) {

  p <- ggplot2::ggplot()
  data <- data.frame()
  groups <- unique(x$counts$group)

  for (i in 2:length(groups)) {
    previous <- x$counts %>%
      dplyr::filter(.data$group == groups[i-1])

    previous_t <- previous$time[nrow(previous)]
    previous_n <- previous$count[nrow(previous)]

    current <- x$counts %>%
      dplyr::filter(.data$group == groups[i])

    if (i == length(groups)) {
      final_t <- if (is.null(final_time)) current$time[nrow(current)] else final_time
    } else {
      final_t <- current$time[nrow(current)]
    }

    xs <- seq(0, final_t - previous_t, length = 100)

    # Extract fit info
    fit_name <- paste0("fit", groups[i])
    fit <- x$fits[[fit_name]]
    ros <- unname(stats::quantile(rstan::extract(fit, pars="ro")$ro, c(.05, .5, .95)))
    K <- mean(rstan::extract(fit, pars="K")$K) * x$fit_info$factor_size

    N_low = K * previous_n / (previous_n + (K - previous_n) * exp(-xs * ros[1]))
    N_medium = K * previous_n / (previous_n + (K - previous_n) * exp(-xs * ros[2]))
    N_high = K * previous_n / (previous_n + (K - previous_n) * exp(-xs * ros[3]))

    d <- data.frame(x=xs + previous_t, yl=N_low, ym=N_medium, yh=N_high)
    p <- p +
      ggplot2::geom_line(d, mapping=ggplot2::aes(x=.data$x, y=.data$ym), col="forestgreen") +
      ggplot2::geom_ribbon(d, mapping=ggplot2::aes(x=.data$x, y=.data$ym, ymin=.data$yl, ymax=.data$yh), fill="forestgreen", alpha=.5)
  }

  p <- p +
    ggplot2::geom_point(x$counts, mapping=ggplot2::aes(x=.data$time, y=.data$count)) +
    my_ggplot_theme()

  if (add_title) {
    p <- p + ggplot2::labs(
      title = paste("Logistic fit", x$sample),
      x = "time",
      y = "count"
    )
  } else {
    p <- p + ggplot2::labs(
      x = "time",
      y = "count"
    )
  }

  p
}

exp_growth = function(t, t0, t_array, rho_array) {

  if (length(t_array) == 0) return(exp(rho_array[1] * (t - t0)))

  if (t <= t_array[1]) {
    return(exp(rho_array[1] * (t - t0)))
  }

  res <- exp(rho_array[1] * (t_array[1] - t0))

  if (length(t_array) >= 2) {
    for (i in 2:length(t_array)) {
      if (t <= t_array[i]) {
        return(res * exp(rho_array[i] * (t - t_array[i-1])))
      } else {
        res <- res * exp(rho_array[i] * (t_array[i] - t_array[i-1]))
      }
    }
  }

  res <- res * exp(rho_array[length(rho_array)] * (t - t_array[length(t_array)]))
  return(res)
}



log_growth = function(t, n0, rho, K) {
  num = K * n0
  den = n0 + (K - n0) * exp(-rho * t)
  return(num/den)
}

log_growth_multiple = function(t, t0, t_array, rho_array, K) {

  current_n0 = 1
  if (length(t_array) == 0) return(log_growth(t - t0, current_n0, rho_array[1], K))

  if (t <= t_array[1]) {
    return(log_growth(t - t0, current_n0, rho_array[1], K))
  }

  dt = t_array[1] - t0
  current_n0 = log_growth(dt, current_n0, rho_array[1], K)
  if (length(t_array) >= 2) {
    for (i in 2:length(t_array)) {
      if (t <= t_array[i]) {
        dt = t - t_array[i-1]
        return(log_growth(dt, current_n0, rho_array[i], K))
      } else {
        dt = t_array[i] - t_array[i-1]
        current_n0 = log_growth(dt, current_n0, rho_array[i], K)
      }
    }
  }

  dt = t - t_array[length(t_array)]
  return(log_growth(dt, current_n0, rho_array[length(rho_array)], K))
}
