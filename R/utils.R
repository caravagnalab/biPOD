
parse_stan_fit = function(f) {
  summary = f$summary() %>% dplyr::as_tibble()
  draws_list = f$draws(format = "draws_list")
  list(summary=summary, draws = draws_list)
}

