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
init <- function(counts, sample, break_points = NULL) {
  cli::cli_h1("biPOD - bayesian inference for Population Dynamics")
  cat("\n")

  # Output
  bipod <- list()
  class(bipod) <- "bipod"

  # Add sample to metadata
  bipod$metadata <- list(sample = sample)
  cli::cli_alert_info("Using sample named: {.field {sample}}.")

  # Parse input
  counts <- check_input_data(counts)
  bipod$counts <- counts

  # Convert breakpoints to groups
  if (is.null(break_points)) {
    bipod$counts$group <- rep(0, nrow(counts))
  } else {
    break_points <- check_break_points(d = counts, break_points = break_points)
    bipod$counts$group <- bp_to_groups(counts, break_points)
  }

  bipod$metadata$breakpoints <- break_points

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

  if (!("group" %in% names(counts))) {
    cli::cli_alert_warning("No group column present in input dataframe! A column will be added.")
    counts <- counts %>%
      dplyr::mutate(group = 0)
  }

  counts <- counts %>%
    dplyr::select(.data$time, .data$count, .data$group) %>%
    dplyr::as_tibble()
  counts
}

bp_to_groups <- function(d, break_points) {
  j <- 1
  groups <- rep(0, nrow(d))
  times <- d$time
  for (i in 1:length(times)) {
    if (j <= length(break_points)) {
      if (times[i] > break_points[j]) j <- j + 1
    }
    groups[i] <- j - 1
  }

  return(groups)
}

check_break_points <- function(d, break_points) {
  if (is.unsorted(break_points)) stop("Break points should be in increasing order")

  # takes breakpoints and d and create group column
  if (any(break_points >= max(d$time))) {
    cli::cli_alert_warning("Some breakpoint are higher than maximum observational time. Extra break points will be removed")
    break_points <- break_points[which(break_points < max(d$time))]
  }

  break_points
}
