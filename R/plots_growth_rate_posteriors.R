#' Produces a list of plots, representing the posteriors, one for each growth rate.
#'
#' @param x a bipod object with a 'fit' field
#'
#' @return A ggplot object containing the posterior density plots of the growth rates and the prior density plot
#' @export
#'
plots_of_growth_rate_posteriors <- function(x) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!("fit" %in% names(x))) stop("Input must contain a 'fits' field")

  # POSTERIOR
  # Obtain list of parameters to plot
  n_groups <- length(unique(x$counts$group))
  par_list <- paste0("rho[", c(1:n_groups), "]")
  colors <- biPOD:::get_group_colors()

  plots <- lapply(1:length(par_list), function(i) {
    par_name <- par_list[i]
    color <- colors[i]

    p <- biPOD:::plot_posterior(x = x, x_fit = x$fit, par_name = par_name, color = color)
    p
  })

  return(plots)
}
