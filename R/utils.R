
# check the column of the input data frame
check_input = function(d) {
  if (!("count" %in% names(d) | "pop.size" %in% names(d))) stop("Input must have a 'count' column or a 'pop.size' column")
  if (!("time" %in% names(d))) stop("Input must have a 'time' column")
  if (!("sample" %in% names(d))) stop("Input must have a 'sample' column")
  if (length(unique(d$sample)) > 1) stop("'sample' must be unique")
}

# clean the data frame by selecting only useful columns
clean_input = function(d) {
  if ("pop.size" %in% names(d)) {
    d <- d %>%
      rename(count = pop.size)
  }
  d <- d %>% select(time, count, sample)
  d
}
