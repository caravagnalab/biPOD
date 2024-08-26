#' Simulate a Stochastic Birth-Death Process
#'
#' Simulates a stochastic birth-death process over multiple time steps. The process models population dynamics where both birth and death rates can vary over time.
#'
#' @param n0 Numeric value indicating the initial population size. Must be a non-negative integer.
#' @param lambda Numeric vector of birth rates.
#'  This can be a single value (constant birth rate) or a vector of length `steps` specifying the birth rate at each time step.
#' @param mu Numeric vector of death rates.
#'  This can be a single value (constant death rate) or a vector of length `steps` specifying the death rate at each time step.
#' @param steps Integer value specifying the number of time steps to simulate. Must be a non-negative integer.
#' @param delta_t Numeric vector of time step sizes. This can be a single value (constant time step) or a vector of length `steps`.
#'  If a single value is provided, it is replicated for each time step.
#'
#' @return A `tibble` with columns:
#' \itemize{
#'   \item \code{time}: The cumulative time at each step.
#'   \item \code{count}: The population size at each time step.
#'   \item \code{group}: A factor indicating the group of contiguous intervals where the birth and death rates are constant.
#' }
#'
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

  # Parameters check
  if (!(n0 >= 0)) stop("n0 must be a positive integer")
  if (!(all(lambda >= 0))) stop("lambda must be positive")
  if (!(all(mu >= 0))) stop("mu must be positive")
  if (!(steps >= 0)) stop("steps must an integer >= 0")
  if (!(all(delta_t >= 0))) stop("delta_t must be positive")

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

  # Produce the groups
  groups <- c(0, group_contiguous(lambda - mu))

  # Return the results as a tibble
  d <- dplyr::tibble(time, count = pop.size, group = groups)
  return(d)
}


sim_single_stochastic_exponential <- function(n0, lambda, mu, delta_t) {
  # Parameters check
  if (!(n0 >= 0)) stop("n0 must be a positive integer")
  if (!(lambda >= 0)) stop("lambda must be positive")
  if (!(mu >= 0)) stop("mu must be positive")
  if (!(delta_t >= 0)) stop("delta_t must be positive")

  # Return n0 if n0 is zero
  if (n0 == 0) {
    return(n0)
  }

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
    event_time <- event_time - log(stats::runif(1, 0, 1)) / (pop_size * w)

    # Update the t variable and counter while the event time is still too small
    while (event_time > t && counter <= 1) {
      counter <- counter + 1
      t <- counter * delta_t
    }

    # Generate a random number to determine if a birth or death event occurs
    if (stats::runif(1, 0, 1) <= p) {
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

#' Simulate Population Growth with Logistic Dynamics
#'
#' Simulates population growth over multiple time steps using a stochastic logistic growth model.
#' This model incorporates both stochastic elements and a carrying capacity to describe population dynamics.
#'
#' @param n0 Numeric value indicating the initial population size. Must be a non-negative integer.
#' @param lambda Numeric vector of birth rates.
#'  This can be a single value (constant birth rate) or a vector of length `steps` specifying the birth rate at each time step.
#' @param mu Numeric vector of death rates.
#'  This can be a single value (constant death rate) or a vector of length `steps` specifying the death rate at each time step.
#' @param K Numeric value representing the carrying capacity of the environment. Must be a non-negative integer.
#' @param steps Integer value specifying the number of time steps to simulate. Must be a non-negative integer.
#' @param delta_t Numeric vector of time step sizes. This can be a single value (constant time step) or a vector of length `steps`.
#'  If a single value is provided, it is replicated for each time step.
#'
#' @return A `tibble` with columns:
#' \itemize{
#'   \item \code{time}: The cumulative time at each step.
#'   \item \code{count}: The population size at each time step.
#'   \item \code{group}: A factor indicating the group of contiguous intervals where the birth and death rates are constant.
#' }
#'
#' @export
sim_stochastic_logistic <- function(n0, lambda, mu, K, steps, delta_t) {
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

  # Parameters check
  if (!(n0 >= 0)) stop("n0 must be a positive integer")
  if (!(all(lambda >= 0))) stop("lambda must be positive")
  if (!(all(mu >= 0))) stop("mu must be positive")
  if (!(K >= 0)) stop("K must be a positive integer")
  if (!(steps >= 0)) stop("steps must an integer >= 0")
  if (!(all(delta_t >= 0))) stop("delta_t must be positive")

  # Return the initial population size if no steps are requested
  if (steps == 0) {
    return(c(n0))
  }

  # Initialize the population size vector
  pop.size <- rep(0, steps + 1)
  pop.size[1] <- n0

  # Simulate the stochastic process at each time step
  for (i in c(2:length(pop.size))) {
    pop.size[i] <- sim_single_stochastic_logistic(pop.size[i - 1], lambda[i - 1], mu[i - 1], K, delta_t[i])
  }

  # Calculate the time, birth rates, and death rates vectors
  time <- cumsum(delta_t)

  # Produce the groups
  groups <- c(0, group_contiguous(lambda - mu))

  # Return the results as a tibble
  d <- dplyr::tibble(time, count = pop.size, group = groups)
  return(d)
}

sim_single_stochastic_logistic <- function(n0, lambda, mu, K, delta_t) {
  # Parameters check
  if (!(n0 >= 0)) stop("n0 must be a positive integer")
  if (!(lambda >= 0)) stop("lambda must be positive")
  if (!(mu >= 0)) stop("mu must be positive")
  if (!(K > 0)) stop("K must be a positive integer")
  if (!(delta_t >= 0)) stop("delta_t must be positive")

  # Calculate the intrinsic growth rate, which is the difference between
  # the per-capita birth rate and death rate.
  ro <- lambda - mu

  # Calculate the density-dependent term, which describes
  # how the growth rate changes as the population size approaches the carrying capacity.
  b2 <- ro / K

  # Set the time of the last event to zero and the current population size to the initial population size.
  event_time <- 0
  pop_size <- n0

  # Initialize a counter and a time variable for the simulation loop.
  counter <- 1
  t <- delta_t

  # Run the simulation loop until the population size reaches zero or until the counter exceeds 1.
  while (TRUE) {
    # Calculate the birth and death rates for the current population size.
    birth <- lambda - b2 * pop_size
    death <- mu

    # Calculate the overall rate of change in the population size.
    rate <- (lambda - b2 * pop_size) * pop_size + (death) * pop_size

    # Update the time of the next event by sampling from an exponential distribution with rate parameter "rate".
    event_time <- event_time - log(stats::runif(1, 0, 1)) / rate

    # Update the time variable for the simulation loop until it exceeds the time step or the counter exceeds 1.
    while (event_time > t && counter <= 1) {
      counter <- counter + 1
      t <- counter * delta_t
    }

    # If the population size has reached zero, exit the simulation loop.
    if (pop_size == 0) break

    # Calculate the probability of a birth event occurring at the next event.
    p <- (lambda - b2 * pop_size) * pop_size / ((lambda - b2 * pop_size) * pop_size + (mu) * pop_size)

    # If the probability is not a number (NA), exit the simulation loop.
    if (is.na(p)) break

    # If a random uniform value is less than or equal to p, increment the population size by 1. Otherwise, decrement it by 1.
    if (stats::runif(1, 0, 1) <= p) {
      pop_size <- pop_size + 1
    } else {
      pop_size <- pop_size - 1
    }

    # If the counter exceeds 1 or the population size reaches zero, exit the simulation loop.
    if (counter > 1 || pop_size == 0) break
  }

  # Return the final population size after the simulation.
  pop_size
}
