
#' Plot the fit over the input data.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param legend_labels Vector containing a name for each unique group in x$counts$group
#' @param legend_title Title for the legend. Default is "group"
#' @param zoom_limits Limits of the x-axis for the zoom of the plot
#' @param sec_axis_breaks Vector containing values in which secondary x axis' breaks will be
#' @param sec_axis_labels Vector containing labels for the secondary x axis
#'
#' @returns A plot of the fit over the input data.
#' @export
plot_fit = function(x,
                    legend_labels = NULL,
                    legend_title = "group",
                    zoom_limits = NULL,
                    sec_axis_breaks = NULL,
                    sec_axis_labels = NULL) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  growth_type <- x$fit_info$growth_type

  if (growth_type == "exponential") {
    fitted_data <- get_exponential_data(x)
  } else {
    fitted_data <- get_logistic_data(x)
  }

  # Plot data
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(x$counts, mapping = ggplot2::aes(x=.data$time, y=.data$count)) + #original points
    ggplot2::geom_line(fitted_data, mapping = ggplot2::aes(x=.data$x, y=.data$y), col="black") +
    ggplot2::geom_ribbon(fitted_data, mapping = ggplot2::aes(x=.data$x, y=.data$y, ymin=.data$ylow, ymax=.data$yhigh), fill="black", alpha=.3) +
    my_ggplot_theme()

  # add highlights
  p <- add_shadow_to_plot(x, base_plot = p)

  # add t0 posterior
  p <- add_t0_posterior(base_plot = p, x = x)

  # change legend
  if (!(is.null(legend_labels))) {
    if (!(length(legend_labels) == length(unique(x$counts$group)))) stop("The number of labels for the legend must be equal to the number of groups")
    p <- p + ggplot2::scale_fill_manual(values = get_group_colors(), labels=legend_labels)
  }
  p <- p + ggplot2::guides(fill=ggplot2::guide_legend(title=legend_title, override.aes = list(alpha = 1)))

  # Add zoom
  zoom_limits <- if (is.null(zoom_limits)) c(min(x$counts$time), max(x$counts$time)) else zoom_limits
  p <- p + ggforce::facet_zoom(xlim = zoom_limits)

  # Add secondary axis
  if (xor(is.null(sec_axis_breaks), is.null(sec_axis_labels))) stop("To plot secondary x axis both sec_axis_breaks and sex_axis_labels are needed")
  if (!(is.null(sec_axis_breaks))) {
    p <- p +
      ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~ ., breaks = sec_axis_breaks, labels = sec_axis_labels)) +
      ggplot2::theme(axis.text.x.top = ggplot2::element_text(angle = 90, vjust = 0.5))
  }

  return(p)
}


#' Plot a simple fit over the input data.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param legend_labels Vector containing a name for each unique group in x$counts$group
#' @param legend_title Title for the legend. Default is "group"
#' @param full_process Boolean, indicating whether to plot the process starting from t0 or not.
#'
#' @returns A plot of the fit over the input data.
#' @export
plot_simple_fit = function(x,
                           full_process = F,
                           legend_labels = NULL,
                           legend_title = "group") {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  growth_type <- x$fit_info$growth_type

  if (growth_type == "exponential") {
    fitted_data <- get_exponential_data(x)
  } else {
    fitted_data <- get_logistic_data(x)
  }

  xmin <- min(x$counts$time)
  xmax <- max(x$counts$time)
  if (!(full_process)) fitted_data <- fitted_data %>% dplyr::filter(.data$x <= xmax) %>% dplyr::filter(.data$x >= xmin)

  # Plot data
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(x$counts, mapping=ggplot2::aes(x=.data$time, y=.data$count)) + #original points
    ggplot2::geom_line(fitted_data, mapping=ggplot2::aes(x=.data$x, y=.data$y), col="black") +
    ggplot2::geom_ribbon(fitted_data, mapping=ggplot2::aes(x=.data$x, y=.data$y, ymin=.data$ylow, ymax=.data$yhigh), fill="black", alpha=.3) +
    my_ggplot_theme()

  # add highlights
  p <- add_shadow_to_plot(x, base_plot = p)

  # add t0 posterior
  if (full_process) p <- add_t0_posterior(base_plot = p, x = x)

  # change legend
  if (!(is.null(legend_labels))) {
    if (!(length(legend_labels) == length(unique(x$counts$group)))) stop("The number of labels for the legend must be equal to the number of groups")
    p <- p + ggplot2::scale_fill_manual(values = get_group_colors(), labels=legend_labels)
  }
  p <- p + ggplot2::guides(fill=ggplot2::guide_legend(title=legend_title, override.aes = list(alpha = 1)))

  # Remove everyithing except process eventually
  if (!(full_process)) p <- p + ggplot2::scale_x_continuous(limits = c(xmin, xmax))

  return(p)
}

get_exponential_data = function(x) {
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
  t0_lower_bound <- x$fit_info$t0_lower_bound

  # Plot model with t0
  if (t0_lower_bound == x$counts$time[1]) {
    t0 <- as.array(t0_lower_bound)
  } else {
    t0 <- rstan::extract(fit, pars=c('t0')) %>% as.list() %>% unlist() %>% as.numeric()
    t0 <- round(t0, 2)
  }

  rho <- rstan::extract(fit, pars=c("rho")) %>% dplyr::as_tibble()
  ro_quantiles <- apply(rho, 2, function(x) stats::quantile(x, c(0.05, 0.5, 0.95))) %>% dplyr::as_tibble() %>% t() %>% dplyr::as_tibble()

  mode_t0 <- stats::median(t0)
  # real_t0 <- (log(factor_size) - ro_quantiles[2]*mode_t0) / (-ro_quantiles[2])
  xs <- seq(mode_t0, max(x$counts$time), length=1000)
  ylow <- lapply(xs, exp_growth, t0=mode_t0, t_array = as.array(t_array), rho_array = ro_quantiles$V1) %>% unlist()
  ymid <- lapply(xs, exp_growth, t0=mode_t0, t_array = as.array(t_array), rho_array = ro_quantiles$V2) %>% unlist()
  yhigh <- lapply(xs, exp_growth, t0=mode_t0, t_array = as.array(t_array), rho_array = ro_quantiles$V3) %>% unlist()

  fitted_data <- dplyr::tibble(
    x = xs,
    y = factor_size * ymid,
    ylow = factor_size * ylow,
    yhigh = factor_size * yhigh
  )

  return(fitted_data)
}

get_logistic_data = function(x) {
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
  t0_lower_bound <- x$fit_info$t0_lower_bound

  # Plot model with t0
  if (t0_lower_bound == x$counts$time[1]) {
    t0 <- as.array(t0_lower_bound)
  } else {
    t0 <- rstan::extract(fit, pars=c('t0')) %>% as.list() %>% unlist() %>% as.numeric()
    t0 <- round(t0, 2)
  }

  K <- mean(rstan::extract(fit, pars=c('K')) %>% as.list() %>% unlist())

  rho <- rstan::extract(fit, pars=c("rho")) %>% dplyr::as_tibble()
  ro_quantiles <- apply(rho, 2, function(x) stats::quantile(x, c(0.05, 0.5, 0.95))) %>% dplyr::as_tibble() %>% t() %>% dplyr::as_tibble()

  mode_t0 <- stats::median(t0)
  xs <- seq(mode_t0, max(x$counts$time), length=1000)

  ylow <- lapply(xs, log_growth_multiple, t0=mode_t0, t_array = as.array(t_array), rho_array = ro_quantiles$V1, K=K) %>% unlist()
  ymid <- lapply(xs, log_growth_multiple, t0=mode_t0, t_array = as.array(t_array), rho_array = ro_quantiles$V2, K=K) %>% unlist()
  yhigh <- lapply(xs, log_growth_multiple, t0=mode_t0, t_array = as.array(t_array), rho_array = ro_quantiles$V3, K=K) %>% unlist()

  fitted_data <- dplyr::tibble(
    x = xs,
    y = factor_size * ymid,
    ylow = factor_size * ylow,
    yhigh = factor_size * yhigh
  )

  return(fitted_data)
}

add_t0_posterior = function(base_plot, x) {

  t0_lower_bound <- x$fit_info$t0_lower_bound

  # Plot model with t0
  if (t0_lower_bound == x$counts$time[1]) {

    base_plot <- base_plot +
      ggplot2::geom_segment(data = data.frame(x = t0_lower_bound, ymin = 0, ymax = max(x$counts$count)), mapping = ggplot2::aes(x=.data$x, xend=.data$x, y=.data$ymin, yend=.data$ymax), col = "darkorange")

  } else {
    values <- rstan::extract(x$fit, pars=c('t0')) %>% as.list() %>% unlist() %>% as.numeric()

    # Add t0 posterior
    df <- get_normalized_density(values, max_value = max(x$counts$count))

    base_plot <- base_plot +
      ggplot2::geom_line(data = df, mapping = ggplot2::aes(x=.data$x, y=.data$y), col = "darkorange") +
      ggplot2::geom_ribbon(data = df, mapping = ggplot2::aes(x=.data$x, ymin=0, ymax=.data$y), fill="darkorange", alpha=.5)
  }
  base_plot
}
