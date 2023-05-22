
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

  counts <- counts %>% dplyr::select(.data$time, .data$count, .data$group) %>% dplyr::as_tibble()
  counts
}

bp_to_groups = function(d, break_points) {

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

check_break_points = function(d, break_points) {
  if (is.unsorted(break_points)) stop("Break points should be in increasing order")

  # takes breakpoints and d and create group column
  if (any(break_points >= max(d$time))) {
    cli::cli_alert_warning("Some breakpoint are higher than maximum observational time. Extra break points will be removed")
    break_points <- break_points[which(break_points < max(d$time))]
  }

  break_points
}
