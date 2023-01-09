
#' Simulate a stochastic birth-death process with logistic grwoth for one time step
#'
#' @param n0 Integer value greater than zero. Starting population.
#' @param lambda Real value greater than zero. Birth rate
#' @param mu Real value greater than zero. Death rate
#' @param K Integer value, carrying capacity.
#' @param delta_t Real value strictly greater than zero. Interval of time.
#' @returns An integer number, indicating the population size at the end of the process.
#' @importFrom stats runif
#' @export
sim_single_stochastic_logistic = function(n0, lambda, mu, K, delta_t) {

  ro <- lambda - mu
  b2 <- ro / K

  event_time <- 0
  pop_size <- n0

  counter <- 1
  t <- delta_t

  while(TRUE) {
    birth <- lambda - b2 * pop_size
    death <- mu
    rate <- (lambda - b2 * pop_size) * pop_size + (death) * pop_size

    event_time <- event_time- log(runif(1,0,1)) / rate
    while (event_time > t && counter <= 1) {
      counter <- counter + 1
      t <- counter * delta_t
    }

    if (pop_size == 0) break

    p <- (lambda - b2 * pop_size) * pop_size / ((lambda - b2 * pop_size) * pop_size + (mu) * pop_size)

    if (is.na(p)) break

    if (runif(1, 0, 1) <= p) {
      pop_size <- pop_size + 1
    } else {
      pop_size <- pop_size - 1
    }

    if (counter > 1 || pop_size == 0) break
  }
  pop_size

}

#' Simulate a stochastic birth-death process for multiple time steps
#'
#' @param n0 Integer value greater than zero. Starting population.
#' @param lambda Vector of real values or single real value. Must be greater than zero. Birth rates.
#' @param mu Vector of real values or single real value. Must be greater than zero. Death rates.
#' #' @param K Integer value, carrying capacity.
#' @param steps Integer value strictly greater than zero. Number of time steps.
#' @param delta_t Vector of real values or single real value. Must be greater than zero. Time intervals.
#' @returns A data frame with the columns: "time", "count", "death.rates", and "birth.rates".
#' @importFrom dplyr tibble
#' @export
sim_stochastic_logistic <- function(n0, lambda, mu, K, steps, delta_t) {
  if (length(mu) == 1) {
    mu <- rep(mu, steps)
  } else {
    stopifnot(steps == length(mu))
  }

  if (length(lambda) == 1) {
    lambda <- rep(lambda, steps)
  } else {
    stopifnot(steps == length(lambda))
  }

  if (length(delta_t) == 1) {
    delta_t <- c(0, rep(delta_t, steps))
  } else {
    stopifnot(steps == length(lambda))
    delta_t <- c(0, delta_t)
  }

  if (steps == 0) {
    return(c(n0))
  }

  pop.size <- rep(0, steps + 1)
  pop.size[1] <- n0

  for (i in c(2:length(pop.size))) {
    pop.size[i] <- sim_single_stochastic_logistic(pop.size[i - 1], lambda[i - 1], mu[i - 1], K, delta_t[i])
  }

  # produce final results
  time <- cumsum(delta_t)
  birth.rates <- c(lambda, 0)
  death.rates <- c(mu, 0)

  d <- tibble(time, pop.size, birth.rates, death.rates)
  return(d)
}
