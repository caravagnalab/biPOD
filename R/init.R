#' Creates a biPOD object.
#'
#' @description Creates a biPOD object from a table containing
#' population counts ad different time points and sample name.
#'
#' @param counts A dataframe of counts with the following fields:
#'
#' * `time` time step;
#' * `count` population count, integer. Alternatively can be named `pop.size`;
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
  assertthat::assert_that("count" %in% names(counts) | "pop.size" %in% names(counts), msg = "input dataframe should contain a column named either 'count' or 'pop.size'")
  assertthat::assert_that("time" %in% names(counts), msg = "input dataframe should contain a column named 'time'")
  assertthat::assert_that(all(counts$count >= 0), msg = "The values of the 'count' column should be all greater or equal than zero!")

  if (is.unsorted(counts$time)) {
    cli::cli_alert_danger("The input should be sorted according to time!")
  }

  if ("pop.size" %in% names(counts)) {
    counts <- counts %>%
      dplyr::rename(count = .data$pop.size)
  }

  if ("group" %in% names(counts)) {
    cli::cli_alert_info("Input sample contains {length(unique(d$group))} group{?s}")
    if (is.unsorted(counts$group)) {
      cli::cli_alert_danger("The groups should be sorted in increasing order!")
    }
  } else {
    cli::cli_alert_warning("Input sample does not specify different groups. A unique group will be considered.")
    counts$group <- rep(0, nrow(counts))
  }
  counts$group[1] <- -1

  counts <- counts %>% dplyr::select(.data$time, .data$count, .data$group)
  counts
}
