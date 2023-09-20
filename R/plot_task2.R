#' Plot the posterior distributions of the inferred breakpoints.
#'
#' @param x a bipod object with a 'fit' field
#' @param with_histogram .
#' @param alpha .
#' @param colors colors to use for the different inferred breakpoints
#'
#' @return A ggplot object containing the posterior distributions of the inferred breakpoints.
#' @export
#'
plot_breakpoints_posterior <- function(x, with_histogram = F, alpha = .6, colors = NULL) {
  x <- x

  # FUN
  # plot_breakpoints_posterior function
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("breakpoints_fit" %in% names(x))) stop("Input must contain a 'breakpoints_fit' field")


  n_changepoints <- length(x$metadata$breakpoints)
  breakpoints_names <- lapply(1:n_changepoints, function(i) {
    paste0("b[", i, "]")
  }) %>% unlist()

  plot_posteriors(x, x$breakpoints_fit,
    par_list = breakpoints_names,
    with_histogram = with_histogram,
    alpha = alpha,
    colors = colors
  )
}
