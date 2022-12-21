
prepare_input_data <- function(counts) {
  stopifnot(("count" %in% names(counts) | "pop.size" %in% names(counts)))
  stopifnot("time" %in% names(counts))

  if ("pop.size" %in% names(counts)) {
    counts <- counts %>%
      dplyr::rename(count = .data$pop.size)
  }
  counts <- counts %>% dplyr::select(.data$time, .data$count)
  counts
}
