#' Create a biPOD Object
#'
#' Initializes a biPOD object from a data frame containing population counts at different time points, along with a sample name and optional breakpoints.
#'
#' @description This function creates a `bipod` object for Bayesian inference of population dynamics.
#'  It takes in a data frame of population counts, a sample name, and optional breakpoints to define groups within the data.
#'
#' @param counts A data frame with two columns:
#' - `time`: Numeric or integer values representing the time steps at which population counts were recorded.
#' - `count`: Integer values representing the population count at each time step.
#'
#' @param sample A character string specifying the name of the sample. This name is stored in the metadata of the resulting biPOD object.
#' @param break_points A numeric vector specifying the breakpoints that define changes in the population dynamics.
#'  If provided, these breakpoints are used to group the time steps. If `NULL`, no grouping is applied. (default is NULL)
#'
#' @return A `bipod` object of class `bipod`, which includes:
#' - A `metadata` list containing the sample name and breakpoints (if provided).
#' - A `counts` data frame with an additional `group` column, which indicates the grouping of time steps based on the provided breakpoints.
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
    min_group_numerosity <- bipod$counts$group %>% table() %>% min()
    # if (min_group_numerosity <= 1) { stop("With the given breakpoints some time windows contain less than 2 observations, which makes the inference not possible") }
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
