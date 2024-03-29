% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_task2.R
\name{breakpoints_inference}
\alias{breakpoints_inference}
\title{Fit growth model to bipod object}
\usage{
breakpoints_inference(
  x,
  factor_size = 1,
  available_changepoints = c(0:3),
  min_support_points = 2,
  max_iter = 20
)
}
\arguments{
\item{x}{a bipod object}

\item{factor_size}{numeric factor by which to divide counts in the bipod object}

\item{available_changepoints}{.}

\item{min_support_points}{.}

\item{max_iter}{.}
}
\value{
the input bipod object with an added 'breakpoints_fit' slot containing the fitted model for the breakpoints
}
\description{
Fit growth model to bipod object
}
