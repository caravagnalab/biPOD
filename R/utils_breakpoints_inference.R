prep_data_bp_inference <- function(x, factor_size = 1, dt = NULL, lambda = .5, n_nodes = NULL) {
  # Get spline prior for changing times
  xs <- x$counts$time
  ys <- x$counts$count

  if (!(is.null(n_nodes))) {
    tuned_lambda <- tune_lambda(xs, ys, n_nodes = n_nodes)
    if (is.null(tuned_lambda)) {
      return(NULL)
    } # not found
    change_x <- spline_nodes(xs, ys, tuned_lambda)
  } else {
    change_x <- spline_nodes(xs, ys, lambda)
  }

  if (is.null(dt)) {
    dt <- min(diff(change_x)) / 2
  }

  # Prepare input data list
  input_data <- list(
    S = nrow(x$counts),
    G = length(as.array(change_x)) + 1,
    N = as.array(as.integer(x$counts$count / factor_size)),
    T = as.array(x$counts$time),
    changing_times_prior = as.array(change_x),
    dt = dt,
    control = list(max_treedepth = 15, adapt_delta = .8)
  )

  return(input_data)
}


tune_lambda <- function(xs, ys, n_nodes) {
  upper <- 1
  lower <- 0

  n_lower <- n_spline_nodes(xs, ys, lower)
  n_upper <- n_spline_nodes(xs, ys, upper)

  if ((n_lower > n_nodes) & (n_upper > n_nodes)) {
    return(NULL)
  }
  if ((n_lower < n_nodes) & (n_upper < n_nodes)) {
    return(NULL)
  }

  while (TRUE) {
    mid <- (lower + upper) / 2 # Calculate the midpoint of the interval
    if (n_spline_nodes(xs, ys, mid) == n_nodes) {
      return(mid)
    } else if (n_spline_nodes(xs, ys, mid) > n_nodes) {
      lower <- mid # Update the lower bound if f(mid) < n
    } else {
      upper <- mid # Update the upper bound otherwise
    }
  }
}

spline_nodes <- function(xs, ys, lambda, n_points = 1000) {
  spl <- stats::smooth.spline(x = xs, y = ys, spar = lambda)

  x_smooth <- seq(min(xs), max(xs), length = n_points)
  y_smooth <- stats::predict(spl, x = x_smooth)[2] %>% unlist()
  dy_smooth <- stats::predict(spl, x = x_smooth, deriv = 1)[2] %>% unlist()

  dy <- diff(y_smooth) / diff(x_smooth)
  dy <- c(dy[1], dy)

  sign_dy <- sign(dy)
  change_idx <- lapply(2:length(dy), function(i) {
    if (sign_dy[i - 1] != sign_dy[i]) {
      return(i)
    }
  }) %>%
    unlist() %>%
    stats::na.omit()

  x_smooth[change_idx]
}

n_spline_nodes <- function(xs, ys, lambda, n_points = 1000) {
  length(spline_nodes(xs, ys, lambda, n_points))
}
