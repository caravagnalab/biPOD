
get_relapse_time_distribution = function(x, n_thresh) {
  # Extract the distribution of the relapse times

  # Extract info from fit
  fit <- x$fit
  G <- length(unique(x$counts$group))

  if (G == 1) {
    t_array = array(0, dim=c(0))
  } else {
    n <- G - 1
    t_array <- (x$counts %>% dplyr::group_by(.data$group) %>% dplyr::slice_tail(n=1) %>% dplyr::ungroup() %>%  dplyr::select(.data$time))$time
    t_array <- t_array[1:n]
  }

  last_t <- max(x$counts$time)

  t0 <- rstan::extract(fit, pars=c('t0')) %>% as.list() %>% unlist() %>% as.numeric()
  t0 <- round(t0, 2)
  best_t0 <- stats::median(t0)
  n_thresh <- n_thresh / x$fit_info$factor_size

  rho <- rstan::extract(fit, pars=c("rho")) %>% dplyr::as_tibble()

  if (x$fit_info$growth_type == "logistic") {
    K <- mean(rstan::extract(fit, pars="K")$K)
    if (K < n_thresh) stop("n_thresh is greater than the predicted carrying capacity")

    t <- max(x$counts$time)
    ys <- lapply(c(1:nrow(rho)), function(i) {
      rho_row <- rho[i,] %>% unlist() %>% unname()

      log_growth_multiple(t = t, t0 = best_t0, t_array = as.array(t_array), rho_array = as.array(rho_row), K = K)
    }) %>% unlist()

    last_rho <- rho[,ncol(rho)] %>% unlist() %>% as.array()

    ys_05 <- stats::quantile(ys, .05)
    ys_95 <- stats::quantile(ys, .95)

    #idx_in_ci <- (ys_05 <= ys) & (ys <= ys_95)

    #ys <- ys[idx_in_ci]
    #last_rho <- last_rho[idx_in_ci]

    #check_last_rho_and_n_thresh(last_rho, ys, n_thresh)
    times <- - (1 / last_rho) * log(ys * (K - n_thresh) / (n_thresh * (K - ys)))

    relapse_times <- times + last_t
  } else {
    # extract the rhos samples

    t <- max(x$counts$time)
    ys <- lapply(1:nrow(rho), function(i) {
      rho_row <- rho[i,] %>% unlist() %>% unname()
      exp_growth(t = t, t0 = best_t0, t_array = as.array(t_array), rho_array = as.array(rho_row))
    }) %>% unlist()

    last_rho <- rho[,ncol(rho)] %>% unlist() %>% as.array()

    ys_05 <- stats::quantile(ys, .05)
    ys_95 <- stats::quantile(ys, .95)

    # idx_in_ci <- (ys_05 <= ys) & (ys <= ys_95)

    # ys <- ys[idx_in_ci]
    #last_rho <- last_rho[idx_in_ci]

    # check_last_rho_and_n_thresh(last_rho = last_rho, last_n = ys, n_thresh = n_thresh)
    relapse_times <- 1 / last_rho * (last_rho * t + log(n_thresh / ys))
  }

  return(dplyr::as_tibble(relapse_times))
}

#' Plot the the posterior distribution of the 'relapse' times.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param n_thresh A positive value indicating the threshold for which
#'  the posterior distribution must be computed
#' @param add_title Boolean, indicating whether the plot should have a title
#'
#' @returns A posterior distribution plot.
#' @export
plot_relapse_time_distribution = function(x, n_thresh, add_title = F) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  times <- get_relapse_time_distribution(x, n_thresh) %>%
    dplyr::mutate(variable = "Relapse_time")

  p <- ggplot2::ggplot(times, ggplot2::aes(x=.data$value)) +
    # ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 100, alpha = .3) +
    ggplot2::geom_density(fill = "violetred4", col = "black", size = .8, alpha = .6) +
    ggplot2::facet_wrap( ~ .data$variable, labeller = ggplot2::label_parsed) +
    ggplot2::labs(
      x = "time",
      y = "density") +
    my_ggplot_theme() +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 4))

  if (add_title) {
    p <- p +
      ggplot2::labs(
        x = "time",
        y = "density",
        title = paste0("Threshold = ", n_thresh)
    )
  } else {
    p <- p +
      ggplot2::labs(
        x = "time",
        y = "density"
      )
  }

  p
}

check_last_rho_and_n_thresh = function(last_rho, last_n, n_thresh) {
  if (all(last_rho <= 0)) {
    # the population is shrinking
    if (!(all(last_n >= n_thresh))) {
      print("RHO INFO")
      print(c(min(last_rho), max(last_rho)))

      print("POPULATION INFO")
      print(c(min(last_n), max(last_n)))
      stop("last inferred growth rate is negative, meaning the population is shrinking, hence the input value of n_thresh will not be reached")
    }
  } else if (all(last_rho >= 0)) {
    # the population is growing
    if (!(all(last_n <= n_thresh))) {
      print("RHO INFO")
      print(c(min(last_rho), max(last_rho)))

      print("POPULATION INFO")
      print(c(min(last_n), max(last_n)))

      stop("last inferred growth rate is positive, meaning the population is expanding, hence the input value of n_thresh will not be reached")
    }
  } else {
    stop("the last inferred rho is ambiguous")
  }
}
