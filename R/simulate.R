#' Simulate a Stochastic Birth-Death Process (Exponential Growth)
#'
#' Simulates a stochastic birth-death process over multiple time steps, which models population dynamics in an environment where both birth and death rates can vary over time.
#' This process follows an exponential growth model, where population size fluctuates due to random birth and death events.
#'
#' @param n0 Numeric value indicating the initial population size at time 0. This must be a non-negative integer (>= 0).
#' @param lambda Numeric vector or scalar specifying the birth rates at each time step.
#' If a scalar is provided, it is assumed to be constant for all time steps. If a vector is provided, it should be of length `steps`, specifying the birth rate at each step.
#' @param mu Numeric vector or scalar specifying the death rates at each time step.
#' Like `lambda`, if a scalar is provided, it is assumed to be constant for all time steps. If a vector is provided, it should be of length `steps`, specifying the death rate at each step.
#' @param steps Integer specifying the number of time steps to simulate. Must be a non-negative integer.
#' @param delta_t Numeric vector or scalar specifying the time step sizes.
#' If a scalar is provided, it is assumed to be constant for all steps. If a vector is provided, it should be of length `steps`, specifying the time step size at each step.
#'
#' @return A `tibble` with the following columns:
#' \itemize{
#'   \item \code{time}: The cumulative time at each step, starting from 0.
#'   \item \code{count}: The population size at each time step, calculated after applying the stochastic birth-death process.
#'   \item \code{group}: A factor that groups contiguous intervals where the birth and death rates are constant (i.e., `lambda` and `mu` remain the same over consecutive time steps).
#' }
#'
#' @details
#' This function simulates a population following an exponential growth model, where the population size changes due to random births and deaths. The rates of birth (`lambda`) and death (`mu`) may vary over time.
#' At each time step, the function computes the expected number of births and deaths based on the rates `lambda` and `mu`. A random draw is then made to determine whether a birth or death event occurs. The population size is updated accordingly.
#'
#' The simulation assumes that the population size is never negative, so if a death event causes the population to drop to zero, it will remain zero thereafter.
#'
#' If `steps` is set to 0, the function returns the initial population size `n0` immediately, without any simulation.
#'
#' @examples
#' # Simulate a stochastic birth-death process with constant rates
#' sim_stochastic_exponential(n0 = 10, lambda = 0.2, mu = 0.1, steps = 10, delta_t = 1)
#'
#' # Simulate with varying rates
#' lambda_rates <- c(rep(0.3, 5), rep(0.1, 5), rep(0.9, 5))
#' mu_rates <- c(rep(0.1, 5), rep(0.15, 5), rep(0.3, 5))
#' sim_stochastic_exponential(n0 = 10, lambda = lambda_rates, mu = mu_rates, steps = 15, delta_t = .1)
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

#' Simulate Population Growth with Stochastic Logistic Dynamics
#'
#' Simulates population growth over multiple time steps using a stochastic logistic growth model.
#' This model incorporates both stochastic elements and a carrying capacity to describe population dynamics.
#' Population size changes due to random births and deaths, but the growth rate is constrained by the carrying capacity \(K\).
#'
#' @param n0 Numeric value indicating the initial population size at time 0. This must be a non-negative integer (>= 0).
#' @param lambda Numeric vector or scalar specifying the birth rates at each time step.
#' If a scalar is provided, it is assumed to be constant for all time steps. If a vector is provided, it should be of length `steps`, specifying the birth rate at each step.
#' @param mu Numeric vector or scalar specifying the death rates at each time step.
#' Like `lambda`, if a scalar is provided, it is assumed to be constant for all time steps. If a vector is provided, it should be of length `steps`, specifying the death rate at each step.
#' @param K Numeric value indicating the carrying capacity of the environment. This must be a non-negative integer (> 0).
#' @param steps Integer specifying the number of time steps to simulate. Must be a non-negative integer.
#' @param delta_t Numeric vector or scalar specifying the time step sizes.
#' If a scalar is provided, it is assumed to be constant for all steps. If a vector is provided, it should be of length `steps`, specifying the time step size at each step.
#'
#' @return A `tibble` with the following columns:
#' \itemize{
#'   \item \code{time}: The cumulative time at each step, starting from 0.
#'   \item \code{count}: The population size at each time step, updated based on the stochastic logistic model.
#'   \item \code{group}: A factor that groups contiguous intervals where the birth and death rates are constant (i.e., `lambda` and `mu` remain the same over consecutive time steps).
#' }
#'
#' @details
#' This function simulates a population following a stochastic logistic growth model, where the population size changes due to random births and deaths while being constrained by the carrying capacity \(K\).
#' The population's growth rate decreases as the population size approaches \(K\), representing the limitations of resources. The rates of birth (`lambda`) and death (`mu`) may vary over time.
#'
#' The stochastic component of the model introduces randomness, making the population size fluctuate between time steps.
#'
#' At each time step, the birth and death rates are used to calculate the expected population change, and a random event (birth or death) is chosen based on the probabilities at that step. The population is updated accordingly.
#'
#' The simulation stops when the population reaches zero or when the specified number of steps is completed.
#'
#' If `steps` is set to 0, the function returns the initial population size `n0` immediately, without any simulation.
#'
#' @examples
#' # Simulate a stochastic logistic growth process with constant rates and carrying capacity of 100
#' sim_stochastic_logistic(n0 = 10, lambda = 0.2, mu = 0.1, K = 100, steps = 15, delta_t = .25)
#'
#' # Simulate with varying rates
#' lambda_rates <- c(rep(0.3, 5), rep(0.1, 5), rep(0.9, 5))
#' mu_rates <- c(rep(0.1, 5), rep(0.15, 5), rep(0.3, 5))
#' sim_stochastic_logistic(
#' n0 = 10,
#' lambda = lambda_rates,
#' mu = mu_rates,
#' steps = 15,
#' delta_t = .1,
#' K=200
#' )
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

group_contiguous <- function(x) {
  rle_x <- rle(x)
  return(rep(seq_along(rle_x$lengths) - 1, times = rle_x$lengths))
}