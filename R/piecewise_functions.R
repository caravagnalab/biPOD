
mean_piecewise_exponential <- function(t, t0, n0, t_array, rho) {
  res <- n0
  if (length(t_array) == 0) return(res * exp(pmax(pmin(rho[1] * (t - t0), 700), -700)))
  if (t <= t_array[1]) return(res * exp(pmax(pmin(rho[1] * (t - t0), 700), -700)))
  res <- res * exp(pmax(pmin(rho[1] * (t_array[1] - t0), 700), -700))
  for (i in 2:length(t_array)) {
    if (t <= t_array[i]) return(res * exp(pmax(pmin(rho[i] * (t - t_array[i - 1]), 700), -700)))
    res <- res * exp(pmax(pmin(rho[i] * (t_array[i] - t_array[i - 1]), 700), -700))
  }
  res * exp(pmax(pmin(rho[length(rho)] * (t - t_array[length(t_array)]), 700), -700))
}

mean_piecewise_logistic <- function(t, t0, n0, t_array, rho, L) {
  res <- n0
  L <- max(L, res * 1.01)  # Ensure L > n0
  if (length(t_array) == 0) {
    exp_term <- exp(pmax(pmin(-rho[1] * (t - t0), 700), -700))
    return(L / (1 + ((L - res) / res) * exp_term))
  }
  if (t <= t_array[1]) {
    exp_term <- exp(pmax(pmin(-rho[1] * (t - t0), 700), -700))
    return(L / (1 + ((L - res) / res) * exp_term))
  }
  exp_term <- exp(pmax(pmin(-rho[1] * (t_array[1] - t0), 700), -700))
  res <- L / (1 + ((L - res) / res) * exp_term)
  for (i in 2:length(t_array)) {
    if (t <= t_array[i]) {
      exp_term <- exp(pmax(pmin(-rho[i] * (t - t_array[i - 1]), 700), -700))
      return(L / (1 + ((L - res) / res) * exp_term))
    }
    exp_term <- exp(pmax(pmin(-rho[i] * (t_array[i] - t_array[i - 1]), 700), -700))
    res <- L / (1 + ((L - res) / res) * exp_term)
  }
  exp_term <- exp(pmax(pmin(-rho[length(rho)] * (t - t_array[length(t_array)]), 700), -700))
  L / (1 + ((L - res) / res) * exp_term)
}

mean_piecewise_gompertz <- function(t, t0, n0, t_array, rho, K) {
  res <- n0
  K <- max(K, res * 1.01)  # Ensure K > n0
  if (length(t_array) == 0) {
    exp_term <- exp(pmax(pmin(-rho[1] * (t - t0), 700), -700))
    return(K * exp(log(pmax(res / K, 1e-10)) * exp_term))
  }
  if (t <= t_array[1]) {
    exp_term <- exp(pmax(pmin(-rho[1] * (t - t0), 700), -700))
    return(K * exp(log(pmax(res / K, 1e-10)) * exp_term))
  }
  exp_term <- exp(pmax(pmin(-rho[1] * (t_array[1] - t0), 700), -700))
  res <- K * exp(log(pmax(res / K, 1e-10)) * exp_term)
  for (i in 2:length(t_array)) {
    if (t <= t_array[i]) {
      exp_term <- exp(pmax(pmin(-rho[i] * (t - t_array[i - 1]), 700), -700))
      return(K * exp(log(pmax(res / K, 1e-10)) * exp_term))
    }
    exp_term <- exp(pmax(pmin(-rho[i] * (t_array[i] - t_array[i - 1]), 700), -700))
    res <- K * exp(log(pmax(res / K, 1e-10)) * exp_term)
  }
  exp_term <- exp(pmax(pmin(-rho[length(rho)] * (t - t_array[length(t_array)]), 700), -700))
  K * exp(log(pmax(res / K, 1e-10)) * exp_term)
}

mean_piecewise_monomolecular <- function(t, t0, n0, t_array, rho, A) {
  res <- n0
  A <- max(A, res * 1.01)  # Ensure A > n0
  if (length(t_array) == 0) {
    exp_term <- exp(pmax(pmin(-rho[1] * (t - t0), 700), -700))
    return(A - (A - res) * exp_term)
  }
  if (t <= t_array[1]) {
    exp_term <- exp(pmax(pmin(-rho[1] * (t - t0), 700), -700))
    return(A - (A - res) * exp_term)
  }
  exp_term <- exp(pmax(pmin(-rho[1] * (t_array[1] - t0), 700), -700))
  res <- A - (A - res) * exp_term
  for (i in 2:length(t_array)) {
    if (t <= t_array[i]) {
      exp_term <- exp(pmax(pmin(-rho[i] * (t - t_array[i - 1]), 700), -700))
      return(A - (A - res) * exp_term)
    }
    exp_term <- exp(pmax(pmin(-rho[i] * (t_array[i] - t_array[i - 1]), 700), -700))
    res <- A - (A - res) * exp_term
  }
  exp_term <- exp(pmax(pmin(-rho[length(rho)] * (t - t_array[length(t_array)]), 700), -700))
  A - (A - res) * exp_term
}

mean_piecewise_quadraticexp <- function(t, t0, n0, t_array, rho) {
  res <- n0
  if (length(t_array) == 0) {
    quad_term <- pmax(pmin(rho[1] * (t - t0)^2, 700), -700)
    return(res * exp(quad_term))
  }
  if (t <= t_array[1]) {
    quad_term <- pmax(pmin(rho[1] * (t - t0)^2, 700), -700)
    return(res * exp(quad_term))
  }
  quad_term <- pmax(pmin(rho[1] * (t_array[1] - t0)^2, 700), -700)
  res <- res * exp(quad_term)
  for (i in 2:length(t_array)) {
    if (t <= t_array[i]) {
      quad_term <- pmax(pmin(rho[i] * (t - t_array[i - 1])^2, 700), -700)
      return(res * exp(quad_term))
    }
    quad_term <- pmax(pmin(rho[i] * (t_array[i] - t_array[i - 1])^2, 700), -700)
    res <- res * exp(quad_term)
  }
  quad_term <- pmax(pmin(rho[length(rho)] * (t - t_array[length(t_array)])^2, 700), -700)
  res * exp(quad_term)
}


# Vectorized piecewise exponential function (corrected)
mean_piecewise_exponential_vec <- function(t, t0, n0, t_array, rho) {
  t <- as.numeric(t)  # Ensure t is numeric vector
  n <- length(t)
  result <- rep(n0, n)

  if (length(t_array) == 0) {
    exp_term <- pmax(pmin(rho[1] * (t - t0), 700), -700)
    return(result * exp(exp_term))
  }

  # Find which segment each t belongs to
  segments <- findInterval(t, c(-Inf, t_array, Inf))

  for (seg in unique(segments)) {
    mask <- segments == seg
    if (!any(mask)) next

    t_subset <- t[mask]

    if (seg == 1) {
      # Before first breakpoint
      exp_term <- pmax(pmin(rho[1] * (t_subset - t0), 700), -700)
      result[mask] <- result[mask] * exp(exp_term)
    } else if (seg <= length(t_array)) {
      # Between breakpoints
      cumulative_growth <- n0

      # Apply growth up to previous breakpoint
      if (seg > 1) {
        exp_term <- pmax(pmin(rho[1] * (t_array[1] - t0), 700), -700)
        cumulative_growth <- cumulative_growth * exp(exp_term)

        if ((seg-1) >= 2) {
          for (i in 2:(seg-1)) {
            exp_term <- pmax(pmin(rho[i] * (t_array[i] - t_array[i-1]), 700), -700)
            cumulative_growth <- cumulative_growth * exp(exp_term)
          }
        }
      }

      # Apply growth in current segment
      prev_t <- if (seg == 1) t0 else t_array[seg-1]
      exp_term <- pmax(pmin(rho[seg] * (t_subset - prev_t), 700), -700)
      result[mask] <- cumulative_growth * exp(exp_term)

    } else {
      # After last breakpoint
      cumulative_growth <- n0

      # Apply all previous growth
      exp_term <- pmax(pmin(rho[1] * (t_array[1] - t0), 700), -700)
      cumulative_growth <- cumulative_growth * exp(exp_term)

      for (i in 2:length(t_array)) {
        exp_term <- pmax(pmin(rho[i] * (t_array[i] - t_array[i-1]), 700), -700)
        cumulative_growth <- cumulative_growth * exp(exp_term)
      }

      # Apply final segment growth
      exp_term <- pmax(pmin(rho[length(rho)] * (t_subset - t_array[length(t_array)]), 700), -700)
      result[mask] <- cumulative_growth * exp(exp_term)
    }
  }

  return(result)
}

# Vectorized piecewise logistic function (corrected)
mean_piecewise_logistic_vec <- function(t, t0, n0, t_array, rho, L) {
  t <- as.numeric(t)
  n <- length(t)
  L <- max(L, n0 * 1.01)
  result <- rep(0, n)

  if (length(t_array) == 0) {
    exp_term <- exp(pmax(pmin(-rho[1] * (t - t0), 700), -700))
    return(L / (1 + ((L - n0) / n0) * exp_term))
  }

  segments <- findInterval(t, c(-Inf, t_array, Inf))

  for (seg in unique(segments)) {
    mask <- segments == seg
    if (!any(mask)) next

    t_subset <- t[mask]

    if (seg == 1) {
      exp_term <- exp(pmax(pmin(-rho[1] * (t_subset - t0), 700), -700))
      result[mask] <- L / (1 + ((L - n0) / n0) * exp_term)
    } else if (seg <= length(t_array)) {
      # Between breakpoints
      cumulative_result <- n0

      # Apply growth up to previous breakpoint
      if (seg > 1) {
        exp_term <- exp(pmax(pmin(-rho[1] * (t_array[1] - t0), 700), -700))
        cumulative_result <- L / (1 + ((L - cumulative_result) / cumulative_result) * exp_term)

        if ((seg-1) >= 2) {
          for (i in 2:(seg-1)) {
            exp_term <- exp(pmax(pmin(-rho[i] * (t_array[i] - t_array[i-1]), 700), -700))
            cumulative_result <- L / (1 + ((L - cumulative_result) / cumulative_result) * exp_term)
          }
        }
      }

      # Apply growth in current segment
      prev_t <- if (seg == 1) t0 else t_array[seg-1]
      exp_term <- exp(pmax(pmin(-rho[seg] * (t_subset - prev_t), 700), -700))
      result[mask] <- L / (1 + ((L - cumulative_result) / cumulative_result) * exp_term)

    } else {
      # After last breakpoint
      cumulative_result <- n0

      # Apply all previous growth
      exp_term <- exp(pmax(pmin(-rho[1] * (t_array[1] - t0), 700), -700))
      cumulative_result <- L / (1 + ((L - cumulative_result) / cumulative_result) * exp_term)

      for (i in 2:length(t_array)) {
        exp_term <- exp(pmax(pmin(-rho[i] * (t_array[i] - t_array[i-1]), 700), -700))
        cumulative_result <- L / (1 + ((L - cumulative_result) / cumulative_result) * exp_term)
      }

      # Apply final segment growth
      exp_term <- exp(pmax(pmin(-rho[length(rho)] * (t_subset - t_array[length(t_array)]), 700), -700))
      result[mask] <- L / (1 + ((L - cumulative_result) / cumulative_result) * exp_term)
    }
  }

  return(result)
}

# Vectorized piecewise Gompertz function (corrected)
mean_piecewise_gompertz_vec <- function(t, t0, n0, t_array, rho, K) {
  t <- as.numeric(t)
  n <- length(t)
  K <- max(K, n0 * 1.01)
  result <- rep(0, n)

  if (length(t_array) == 0) {
    exp_term <- exp(pmax(pmin(-rho[1] * (t - t0), 700), -700))
    return(K * exp(log(pmax(n0 / K, 1e-10)) * exp_term))
  }

  segments <- findInterval(t, c(-Inf, t_array, Inf))

  for (seg in unique(segments)) {
    mask <- segments == seg
    if (!any(mask)) next

    t_subset <- t[mask]

    if (seg == 1) {
      exp_term <- exp(pmax(pmin(-rho[1] * (t_subset - t0), 700), -700))
      result[mask] <- K * exp(log(pmax(n0 / K, 1e-10)) * exp_term)
    } else if (seg <= length(t_array)) {
      # Between breakpoints
      cumulative_result <- n0

      # Apply growth up to previous breakpoint
      if (seg > 1) {
        exp_term <- exp(pmax(pmin(-rho[1] * (t_array[1] - t0), 700), -700))
        cumulative_result <- K * exp(log(pmax(cumulative_result / K, 1e-10)) * exp_term)

        if ((seg-1) >= 2) {
          for (i in 2:(seg-1)) {
            exp_term <- exp(pmax(pmin(-rho[i] * (t_array[i] - t_array[i-1]), 700), -700))
            cumulative_result <- K * exp(log(pmax(cumulative_result / K, 1e-10)) * exp_term)
          }
        }
      }

      # Apply growth in current segment
      prev_t <- if (seg == 1) t0 else t_array[seg-1]
      exp_term <- exp(pmax(pmin(-rho[seg] * (t_subset - prev_t), 700), -700))
      result[mask] <- K * exp(log(pmax(cumulative_result / K, 1e-10)) * exp_term)

    } else {
      # After last breakpoint
      cumulative_result <- n0

      # Apply all previous growth
      exp_term <- exp(pmax(pmin(-rho[1] * (t_array[1] - t0), 700), -700))
      cumulative_result <- K * exp(log(pmax(cumulative_result / K, 1e-10)) * exp_term)

      for (i in 2:length(t_array)) {
        exp_term <- exp(pmax(pmin(-rho[i] * (t_array[i] - t_array[i-1]), 700), -700))
        cumulative_result <- K * exp(log(pmax(cumulative_result / K, 1e-10)) * exp_term)
      }

      # Apply final segment growth
      exp_term <- exp(pmax(pmin(-rho[length(rho)] * (t_subset - t_array[length(t_array)]), 700), -700))
      result[mask] <- K * exp(log(pmax(cumulative_result / K, 1e-10)) * exp_term)
    }
  }

  return(result)
}

# Vectorized piecewise monomolecular function (corrected)
mean_piecewise_monomolecular_vec <- function(t, t0, n0, t_array, rho, A) {
  t <- as.numeric(t)
  n <- length(t)
  A <- max(A, n0 * 1.01)
  result <- rep(0, n)

  if (length(t_array) == 0) {
    exp_term <- exp(pmax(pmin(-rho[1] * (t - t0), 700), -700))
    return(A - (A - n0) * exp_term)
  }

  segments <- findInterval(t, c(-Inf, t_array, Inf))

  for (seg in unique(segments)) {
    mask <- segments == seg
    if (!any(mask)) next

    t_subset <- t[mask]

    if (seg == 1) {
      exp_term <- exp(pmax(pmin(-rho[1] * (t_subset - t0), 700), -700))
      result[mask] <- A - (A - n0) * exp_term
    } else if (seg <= length(t_array)) {
      # Between breakpoints
      cumulative_result <- n0

      # Apply growth up to previous breakpoint
      if (seg > 1) {
        exp_term <- exp(pmax(pmin(-rho[1] * (t_array[1] - t0), 700), -700))
        cumulative_result <- A - (A - cumulative_result) * exp_term

        if ((seg-1) >= 2) {
          for (i in 2:(seg-1)) {
            exp_term <- exp(pmax(pmin(-rho[i] * (t_array[i] - t_array[i-1]), 700), -700))
            cumulative_result <- A - (A - cumulative_result) * exp_term
          }
        }
      }

      # Apply growth in current segment
      prev_t <- if (seg == 1) t0 else t_array[seg-1]
      exp_term <- exp(pmax(pmin(-rho[seg] * (t_subset - prev_t), 700), -700))
      result[mask] <- A - (A - cumulative_result) * exp_term

    } else {
      # After last breakpoint
      cumulative_result <- n0

      # Apply all previous growth
      exp_term <- exp(pmax(pmin(-rho[1] * (t_array[1] - t0), 700), -700))
      cumulative_result <- A - (A - cumulative_result) * exp_term

      for (i in 2:length(t_array)) {
        exp_term <- exp(pmax(pmin(-rho[i] * (t_array[i] - t_array[i-1]), 700), -700))
        cumulative_result <- A - (A - cumulative_result) * exp_term
      }

      # Apply final segment growth
      exp_term <- exp(pmax(pmin(-rho[length(rho)] * (t_subset - t_array[length(t_array)]), 700), -700))
      result[mask] <- A - (A - cumulative_result) * exp_term
    }
  }

  return(result)
}

# Vectorized piecewise quadratic exponential function (corrected)
mean_piecewise_quadraticexp_vec <- function(t, t0, n0, t_array, rho) {
  t <- as.numeric(t)
  n <- length(t)
  result <- rep(n0, n)

  if (length(t_array) == 0) {
    quad_term <- pmax(pmin(rho[1] * (t - t0)^2, 700), -700)
    return(result * exp(quad_term))
  }

  segments <- findInterval(t, c(-Inf, t_array, Inf))

  for (seg in unique(segments)) {
    mask <- segments == seg
    if (!any(mask)) next

    t_subset <- t[mask]

    if (seg == 1) {
      quad_term <- pmax(pmin(rho[1] * (t_subset - t0)^2, 700), -700)
      result[mask] <- result[mask] * exp(quad_term)
    } else if (seg <= length(t_array)) {
      # Between breakpoints
      cumulative_growth <- n0

      # Apply growth up to previous breakpoint
      if (seg > 1) {
        quad_term <- pmax(pmin(rho[1] * (t_array[1] - t0)^2, 700), -700)
        cumulative_growth <- cumulative_growth * exp(quad_term)

        if ((seg-1) >= 2) {
          for (i in 2:(seg-1)) {
            quad_term <- pmax(pmin(rho[i] * (t_array[i] - t_array[i-1])^2, 700), -700)
            cumulative_growth <- cumulative_growth * exp(quad_term)
          }
        }
      }

      # Apply growth in current segment
      prev_t <- if (seg == 1) t0 else t_array[seg-1]
      quad_term <- pmax(pmin(rho[seg] * (t_subset - prev_t)^2, 700), -700)
      result[mask] <- cumulative_growth * exp(quad_term)

    } else {
      # After last breakpoint
      cumulative_growth <- n0

      # Apply all previous growth
      quad_term <- pmax(pmin(rho[1] * (t_array[1] - t0)^2, 700), -700)
      cumulative_growth <- cumulative_growth * exp(quad_term)

      for (i in 2:length(t_array)) {
        quad_term <- pmax(pmin(rho[i] * (t_array[i] - t_array[i-1])^2, 700), -700)
        cumulative_growth <- cumulative_growth * exp(quad_term)
      }

      # Apply final segment growth
      quad_term <- pmax(pmin(rho[length(rho)] * (t_subset - t_array[length(t_array)])^2, 700), -700)
      result[mask] <- cumulative_growth * exp(quad_term)
    }
  }

  return(result)
}
