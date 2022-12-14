#' Creates a biPOD object.
#'
#' @description Creates a biPOD object from a table containing
#' population counts ad different time points and sample name.
#'
#' @param counts A dataframe of counts with the following fields:
#'
#' * `time` time step;
#' * `count` population count, integer. Alternatively can be named 'pop.size';
#'
#' @param sample A string containing the sample name.
#'
#' @return A biPOD object of class `bipod`, with S3 methods for printing,
#' plotting and analyzing data.
#'
#' @export
init = function(counts, sample) {
  cli::cli_h1("biPOD - bayesian inference for Population Dynamics")
  cat("\n")

  # Output
  bipod = list()
  class(bipod) <- "bipod"

  # Sample
  bipod$sample <- sample
  cli::cli_alert_info("Using sample named: {.field {sample}}.")

  # Parse input
  input <- prepare_input_data(counts)
  bipod$counts <- input

  bipod
}
