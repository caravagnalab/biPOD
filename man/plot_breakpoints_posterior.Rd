% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_task0.R
\name{plot_breakpoints_posterior}
\alias{plot_breakpoints_posterior}
\title{Plot Posterior Distributions of Inferred Breakpoints}
\usage{
plot_breakpoints_posterior(x, with_histogram = F, alpha = 0.6, colors = NULL)
}
\arguments{
\item{x}{A \code{bipod} object that contains a 'breakpoints_fit' field.}

\item{with_histogram}{A logical value indicating whether to overlay histograms of the breakpoint samples on the density plots. (default is FALSE)}

\item{alpha}{A numeric value between 0 and 1 specifying the transparency level of the density plot. (default is 0.6)}

\item{colors}{A character vector specifying colors to use for the different inferred breakpoints. If \code{NULL}, the default color 'darkgray' is used for all breakpoints.}
}
\value{
A \code{ggplot2} object displaying the posterior distributions of the inferred breakpoints. The plot shows density estimates and optionally histograms for the breakpoints.
}
\description{
Generates a plot showing the posterior distributions of the inferred breakpoints from a fitted model.
The plot can include density estimates and optional histograms for better visualization.
}
