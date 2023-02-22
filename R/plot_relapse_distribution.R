
get_relapse_time_distribution_old = function(x, n_thresh) {
  # Extract the distribution of the relapse times
  # Extract the last fitted growth rate and the
  # previous count as the n0

  max_group <- max(x$counts$group)

  # Extract 'previous' n0
  previous_counts <- x$counts %>%
    dplyr::filter(.data$group == max_group - 1)
  n0 <- previous_counts$count[nrow(previous_counts)]

  # Extract samples of the last growth rate
  fit_name <- paste0("fit", max_group)
  ro_samples <- rstan::extract(x$fits[[fit_name]], pars=c('ro'))$ro


  # Compute the relapse times
  if (x$fit_info$growth_type == "logistic") {
    K <- mean(rstan::extract(x$fits[[fit_name]], pars="K")$K)
    if (K < n_thresh) stop("n_thresh is greater than the predicted carrying capacity")
    times <- - 1 / ro_samples * log(n0 * (K - n_thresh) / (n_thresh * (K - n0)))
  } else {
    times <- 1 / ro_samples * log(n_thresh / n0)
  }
  return(as.data.frame(times))
}

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
  median_t0 <- median(t0)

  rho <- rstan::extract(fit, pars=c("rho")) %>% as.data.frame()

  if (x$fit_info$growth_type == "logistic") {
    K <- mean(rstan::extract(x$fits[[fit_name]], pars="K")$K)
    if (K < n_thresh) stop("n_thresh is greater than the predicted carrying capacity")

    last_n <- lapply(c(1:nrow(rho)), function(i) {
      rho_array <- rho[i,] %>% as.numeric()
      return(log_growth_multiple(last_t, median_t0, t_array, rho_array, K))
    }) %>% unlist()

    last_rho_idx <- ncol(rho)
    last_rho <- rho[,last_rho_idx]

    check_last_rho_and_n_thresh(last_rho, last_n, n_thresh)
    times <- - (1 / last_rho) * log(last_n * (K - n_thresh) / (n_thresh * (K - last_n)))

  } else {
    # extract the rhos samples

    last_n <- lapply(c(1:nrow(rho)), function(i) {
      rho_array <- rho[i,] %>% as.numeric()
      return(exp_growth(last_t, median_t0, t_array, rho_array))
    }) %>% unlist()

    last_rho_idx <- ncol(rho)
    last_rho <- rho[,last_rho_idx]

    check_last_rho_and_n_thresh(last_rho, last_n, n_thresh)
    times <- (1 / last_rho) * log(n_thresh / last_n)
  }

  times <- times + last_t
  return(as.data.frame(times))
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

  times <- get_relapse_time_distribution(x, n_thresh)
  bw = (max(times$times) - min(times$times)) / 100
  p <- ggplot2::ggplot(times, aes(x=.data$times)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = bw, alpha = .3) +
    ggplot2::geom_density(col = "darkorange", size = .8) +
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

plot_relapse_time_distribution_old = function(x, n_thresh, add_title = F) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  times <- get_relapse_time_distribution(x, n_thresh)
  bw = (max(times$times) - min(times$times)) / 100
  p <- ggplot2::ggplot(times, aes(x=.data$times)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = bw, alpha = .3) +
    ggplot2::geom_density(col = "darkorange", size = .8) +
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
    if (!(all(last_n >= n_thresh))) stop("last inferred growth rate is negative, meaning the population is shrinking, hence the input value of n_thresh will not be reached")
  } else if (all(last_rho >= 0)) {
    # the population is growing
    if (!(all(last_n <= n_thresh))) stop("last inferred growth rate is positive, meaning the population is expanding, hence the input value of n_thresh will not be reached")
  } else {
    stop("the last inferred rho is ambiguous")
  }
}
