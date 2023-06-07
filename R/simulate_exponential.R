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
