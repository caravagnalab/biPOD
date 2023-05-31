
prep_data_bp_inference = function(x, factor_size = 1, prior_K = NULL, dt = NULL, spline_spar=.5) {

  # Parameters check
  if (is.null(prior_K)) {
    prior_K = max(x$counts$count) / factor_size
  } else {
    prior_K = prior_K / factor_size
    if (prior_K <= 0) stop("'prior_K' should eiter be NULL or positive")
  }

  # Get spline prior for changing times
  spl <- stats::smooth.spline(x = x$counts$time, y = x$counts$count, spar = spline_spar)

  x_smooth <- seq(min(x$counts$time), max(x$counts$time), length = 1000)
  y_smooth <- stats::predict(spl, x = x_smooth)[2] %>% unlist()
  dy_smooth <- stats::predict(spl, x = x_smooth, deriv = 1)[2] %>% unlist()

  dy <- diff(y_smooth) / diff(x_smooth)
  dy <- c(dy[1], dy)

  sign_dy <- sign(dy)
  change_idx <- lapply(2:length(dy), function(i) {
    if (sign_dy[i-1] != sign_dy[i]) {return(i)}
  }) %>% unlist() %>% stats::na.omit()

  change_x <- x_smooth[change_idx]

  if (is.null(dt)) { dt <- min(diff(change_x)) / 2}

  # Prepare input data list
  input_data <- list(
    S = nrow(x$counts),
    G = length(as.array(change_x)) + 1,
    N = as.array(as.integer(x$counts$count / factor_size)),
    T = as.array(x$counts$time),
    changing_times_prior = as.array(change_x),
    dt = dt,
    prior_K = prior_K,
    control = list(max_treedepth = 15, adapt_delta = .8)
  )

  return(input_data)
}


