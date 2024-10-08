% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_task3.R
\name{fit_two_pop_model}
\alias{fit_two_pop_model}
\title{Fit a Two-Population Growth Model to a bipod Object}
\usage{
fit_two_pop_model(
  x,
  variational = FALSE,
  factor_size = 1,
  chains = 4,
  iter = 5000,
  cores = 4
)
}
\arguments{
\item{x}{A \code{bipod} object.}

\item{variational}{A logical value indicating whether to use Variational Inference instead of MCMC sampling.
If \code{TRUE}, the model will be fitted using Variational Inference; otherwise, MCMC sampling will be used. (default is FALSE)}

\item{factor_size}{A numeric value representing the factor by which to divide the counts in the bipod object.
This value must be positive and appropriate for the data scale. (default is 1)}

\item{chains}{An integer specifying the number of chains to run in the MCMC algorithm.
This parameter is ignored if \code{variational = TRUE}. (default is 4)}

\item{iter}{An integer specifying the number of iterations to run in the MCMC algorithm.
This parameter is ignored if \code{variational = TRUE}. (default is 5000)}

\item{cores}{An integer specifying the number of cores to use for parallel processing during model fitting. (default is 4)}
}
\value{
The input \code{bipod} object with added slots:
\itemize{
\item \code{'two_pop_fit'}: Contains the fitted two-population growth model.
\item \code{'two_pop_fit_info'}: Contains information about the fitting process, including metadata such as sampling type and factor size.
}
}
\description{
This function fits a growth model that accounts for two populations within a bipod object. The fitting can be done using either Variational Inference or Markov Chain Monte Carlo (MCMC) sampling.
}
