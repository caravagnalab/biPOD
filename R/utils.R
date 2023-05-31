
group_contiguous <- function(x) {
  rle_x <- rle(x)
  return(rep(seq_along(rle_x$lengths) - 1, times = rle_x$lengths))
}

exp_growth = function(t, t0, t_array, rho_array, n0) {

  if (length(t_array) == 0) return(n0 * exp(rho_array[1] * (t - t0)))

  if (t <= t_array[1]) {
    return(n0 * exp(rho_array[1] * (t - t0)))
  }

  res <- n0 * exp(rho_array[1] * (t_array[1] - t0))

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

log_growth_multiple = function(t, t0, t_array, rho_array, K, n0) {

  current_n0 = n0
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

extract_parameter = function(fit, par_name) {
  rstan::extract(fit, pars=c(par_name)) %>% dplyr::as_tibble() %>%
    dplyr::mutate(par_name = par_name) %>%
    dplyr::rename("value" = 1, "parameter" = 2)
}

extract_parameters = function(fit, par_list) {
  d <- dplyr::tibble()
  for (p in par_list) {
    d <- dplyr::bind_rows(d, extract_parameter(fit, p))
  }
  d
}
