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
#'
#' @return A biPOD object of class `bipod`.
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
  input <- check_input_data(counts)
  bipod$counts <- input

  bipod
}

check_input_data <- function(counts) {
  if (!("count" %in% names(counts))) stop("Input dataframe should contain a column named either 'count'")
  if (!("time" %in% names(counts))) stop("Input dataframe should contain a column named either 'time'")
  if (!(all(counts$count >= 0))) stop("The values of the 'count' column should be all be positive or equal to zero")

  if (is.unsorted(counts$time)) {
    cli::cli_alert_danger("The input should be sorted according to time!")
    cli::cli_alert_danger("The input is being sorted according to time...")
    counts <- counts %>%
      dplyr::arrange(.data$time)
  }

  if ("group" %in% names(counts)) {
    cli::cli_alert_info("Input sample contains {length(unique(counts$group))} group{?s}")
    cli::cli_alert_info("The groups are being trasnformed into integer values...")
    counts$group <- group_contiguous(counts$group)
  } else {
    cli::cli_alert_warning("Input sample does not specify different groups. A unique group will be considered.")
    counts$group <- rep(0, nrow(counts))
  }
  cat("\n")

  counts <- counts %>% dplyr::select(.data$time, .data$count, .data$group)
  counts
}
