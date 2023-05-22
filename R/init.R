#' Creates a biPOD object.
#'
#' @description Creates a biPOD object from a table containing
#' population counts ad different time points and sample name.
#'
#' @param counts A dataframe of counts with the following fields:
#'
#' * `time` time step.
#' * `count` population count, integer.
#'
#' @param sample A string containing the sample name.
#' @param break_points An array containing the changing points.
#'
#' @return A biPOD object of class `bipod`.
#'
#' @export
init = function(counts, sample, break_points = NULL) {
  cli::cli_h1("biPOD - bayesian inference for Population Dynamics")
  cat("\n")

  # Output
  bipod = list()
  class(bipod) <- "bipod"

  # Add sample to metadata
  bipod$metadata <- list(sample = sample)
  cli::cli_alert_info("Using sample named: {.field {sample}}.")

  # Parse input
  counts <- biPOD:::check_input_data(counts)
  bipod$counts <- counts

  # Convert breakpoints to groups
  if (is.null(break_points)) {
    bipod$counts$group <- rep(0, nrow(counts))
  } else {
    break_points <- biPOD:::check_break_points(d = counts, break_points = break_points)
    bipod$counts$group <- biPOD:::bp_to_groups(counts, break_points)
  }

  bipod$metadata$breakpoints <- break_points

  bipod
}




