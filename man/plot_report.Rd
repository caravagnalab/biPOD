% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_report.R
\name{plot_report}
\alias{plot_report}
\title{Generate a Comprehensive Report Plot}
\usage{
plot_report(
  x,
  fit_type = "complex",
  breakpoints_color = "darkgray",
  shadows_colors = NULL,
  t0_posterior_color = "darkorange",
  full_process = FALSE
)
}
\arguments{
\item{x}{A \code{bipod} object. Must contain a 'fit' and metadata fields for model selection and, optionally, breakpoints and ELBO values.}

\item{fit_type}{A character string specifying the type of fit plot to generate. Options include "simple" and "complex". (default is "complex")}

\item{breakpoints_color}{A character string specifying the color for the breakpoints posterior plot. (default is "darkgray")}

\item{shadows_colors}{A character vector specifying colors for shadows in the fit plot. (default is NULL)}

\item{t0_posterior_color}{A character string specifying the color for the t0 posterior plot. (default is "darkorange")}

\item{full_process}{A logical value indicating whether to include the full process in the fit plot when \code{fit_type} is "simple". (default is FALSE)}
}
\value{
A \code{patchwork} object containing multiple ggplot2 plots arranged in a single panel.
}
\description{
Creates a detailed report consisting of various plots to visualize the results of a fitted model within a \code{bipod} object.
The report includes fit plots, breakpoints posterior plots, growth rates posterior plots, and model selection plots.
}
