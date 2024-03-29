% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{sim_stochastic_logistic}
\alias{sim_stochastic_logistic}
\title{Simulate population growth with logistic growth for multiple time steps}
\usage{
sim_stochastic_logistic(n0, lambda, mu, K, steps, delta_t)
}
\arguments{
\item{n0}{The initial population size.}

\item{lambda}{The rate at which births occur. Can be a single value or a vector
of length \code{steps}.}

\item{mu}{The rate at which deaths occur. Can be a single value or a vector
of length \code{steps}.}

\item{K}{Carrying capacity.}

\item{steps}{The number of time steps to simulate.}

\item{delta_t}{The size of the time step. Can be a single value or a vector
of length \code{steps}.}
}
\value{
A tibble with the time, population size.
}
\description{
Simulate population growth with logistic growth for multiple time steps
}
