
#' Simulate a single realization of a stochastic process with birth and death events.
#'
#' This function simulates a single realization of a stochastic process where the
#' population size is subject to birth and death events. The rates at which these
#' events occur are given by the `lambda` and `mu` parameters, respectively. The
#' function returns the final population size after the simulation is completed.
#'
#' @param n0 The initial population size.
#' @param lambda The rate at which births occur.
#' @param mu The rate at which deaths occur.
#' @param delta_t The size of the time step.
#' @returns The final population size after the simulation.
#'
#' @importFrom stats runif
#' @export
sim_single_stochastic_exponential <- function(n0, lambda, mu, delta_t) {
  assertthat::assert_that(n0 >= 0 , msg = "n0 must be a positive integer")
  assertthat::assert_that(lambda >= 0 , msg = "lambda must be positive")
  assertthat::assert_that(mu >= 0 , msg = "mu must be positive")
  assertthat::assert_that(delta_t >= 0 , msg = "delta_t must be positive")

  # Calculate the total event rate and the probability of a birth event
  w <- lambda + mu
  p <- lambda / w

  # Initialize the event time and population size
  event_time <- 0
  pop_size <- n0

  # Initialize the counter and t variables for tracking the number of iterations
  counter <- 1
  t <- delta_t

  # Enter an infinite loop
  while (TRUE) {
    # Update the event time by the time until the next event occurs
    event_time <- event_time - log(runif(1, 0, 1)) / (pop_size * w)

    # Update the t variable and counter while the event time is still too small
    while (event_time > t && counter <= 1) {
      counter <- counter + 1
      t <- counter * delta_t
    }

    # Generate a random number to determine if a birth or death event occurs
    if (runif(1, 0, 1) <= p) {
      pop_size <- pop_size + 1
    } else {
      pop_size <- pop_size - 1
    }

    # Break the loop if the counter exceeds 1 or the population size is 0
    if (counter > 1 || pop_size == 0) break
  }

  # Return the final population size
  pop_size
}

#' Simulate a stochastic birth-death process for multiple time steps
#'
#' @param n0 The initial population size.
#' @param lambda The rate at which births occur. Can be a single value or a vector
#'   of length `steps`.
#' @param mu The rate at which deaths occur. Can be a single value or a vector
#'   of length `steps`.
#' @param steps The number of time steps to simulate.
#' @param delta_t The size of the time step. Can be a single value or a vector
#'   of length `steps`.
#' @return A tibble with the time, population size, birth rates, and death rates
#'   at each time step.
#'
#' @importFrom dplyr tibble
#' @export
sim_stochastic_exponential <- function(n0, lambda, mu, steps, delta_t) {
  # Convert single values of lambda, mu, and delta_t to vectors of length steps
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

  # Check input
  assertthat::assert_that(n0 >= 0 , msg = "n0 must be a positive integer")
  assertthat::assert_that(all(lambda >= 0) , msg = "lambda must be positive")
  assertthat::assert_that(all(mu >= 0) , msg = "mu must be positive")
  assertthat::assert_that(steps >= 0 , msg = "steps must a integere greater or equal than 0")
  assertthat::assert_that(all(delta_t >= 0) , msg = "delta_t must be positive")

  # Return the initial population size if no steps are requested
  if (steps == 0) {
    return(c(n0))
  }

  # Initialize the population size vector
  pop.size <- rep(0, steps + 1)
  pop.size[1] <- n0

  # Simulate the stochastic process at each time step
  for (i in c(2:length(pop.size))) {
    pop.size[i] <- sim_single_stochastic_exponential(n0 = pop.size[i - 1], lambda = lambda[i - 1], mu = mu[i - 1], delta_t = delta_t[i])
  }

  # Calculate the time, birth rates, and death rates vectors
  time <- cumsum(delta_t)
  birth.rates <- c(lambda, 0)
  death.rates <- c(mu, 0)

  # Return the results as a tibble
  d <- tibble(time, pop.size, birth.rates, death.rates)
  return(d)
}

#' Simulate population growth with noise
#'
#' This function simulates population growth over a given number of steps, with
#' noise introduced at each step. The population size at each step is calculated
#' based on the growth rate (lambda) and death rate (mu) at that step, and the
#' delta time (delta_t) between steps. Noise is introduced to the population size
#' at each step according to a normal distribution with standard deviation equal
#' to the population size multiplied by a given noise parameter (sigma).
#'
#' @param n0 Initial population size.
#' @param lambda Growth rate. If a single value is provided, it is assumed to be
#'   constant across all steps. If a vector is provided, it must be the same
#'   length as the number of steps.
#' @param mu Death rate. If a single value is provided, it is assumed to be
#'   constant across all steps. If a vector is provided, it must be the same
#'   length as the number of steps.
#' @param steps Number of steps to simulate.
#' @param delta_t Delta time between steps. If a single value is provided, it is
#'   assumed to be constant across all steps. If a vector is provided, it must be
#'   the same length as the number of steps.
#' @param sigma Noise parameter.
#' @param as_int Logical indicating whether to round the population size to the
#'   nearest integer at each step.
#' @return A tibble containing columns for time, population size, birth rate, and
#'   death rate at each step.
#'
#' @importFrom dplyr tibble
#' @importFrom stats rnorm
#' @export
sim_noisy <- function(n0, lambda, mu, steps, delta_t, sigma, as_int = T) {
  # Convert single values of lambda, mu, and delta_t to vectors of length steps
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

  # Check input
  assertthat::assert_that(n0 >= 0 , msg = "n0 must be a positive integer")
  assertthat::assert_that(all(lambda >= 0) , msg = "lambda must be positive")
  assertthat::assert_that(all(mu >= 0) , msg = "mu must be positive")
  assertthat::assert_that(steps >= 0 , msg = "steps must a integere greater or equal than 0")
  assertthat::assert_that(all(delta_t >= 0) , msg = "delta_t must be positive")

  # Return the initial population size if no steps are requested
  if (steps == 0) {
    return(c(n0))
  }

  # Initialize the population size vector
  pop.size <- rep(0, steps + 1)
  pop.size[1] <- n0

  # Simulate the stochastic process at each time step
  for (i in c(2:length(pop.size))) {
    pop.size[i] <- pop.size[i - 1] * exp(delta_t[i] * (lambda[i - 1] - mu[i - 1]))
  }

  # Add noise
  stds <- pop.size * sigma
  noise <- rnorm(length(pop.size), 0, stds)
  pop.size <- pop.size + noise
  if (as_int) {
    pop.size <- as.integer(pop.size)
  }

  # produce final results
  time <- cumsum(delta_t)
  birth.rates <- c(lambda, 0)
  death.rates <- c(mu, 0)
  d <- tibble(time, pop.size, birth.rates, death.rates)

  return(d)
}
