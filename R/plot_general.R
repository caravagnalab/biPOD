#' Plot posterior predictive ribbon
#'
#' Creates a ribbon plot showing the median and credible intervals of the posterior
#' predictive distribution (`y_rep`) extracted from a fitted model.
#'
#' @param fit A list containing posterior draws with element `draws` (list of chains).
#' Optionally, `fit$time` may provide the time points.
#' @param data Optional data frame with columns `time` and `count` to overlay observed data.
#' @param ci Credible interval width (default `0.9`).
#'
#' @return A `ggplot` object showing the posterior predictive ribbon and observed data.
#' @examples
#' \dontrun{
#'   plot_ribbon(fit_object, data = my_data)
#' }
#'
#' @export
plot_ribbon <- function(fit, data = NULL, ci = 0.9, shadow_breakpoints = NULL, shadow_colors = NULL) {

  draws_list <- fit$draws
  if (!is.list(draws_list)) stop("fit$draws must be a list of chains.")

  extract_chain <- function(chain) {
    param_names <- names(chain)
    yrep_pattern <- grep("y[_]?rep\\[", param_names, value = TRUE)
    if (length(yrep_pattern) == 0) stop("No y_rep[...] parameters found.")
    y_rep_mat <- sapply(yrep_pattern, function(p) unlist(chain[[p]]))
    if (is.null(dim(y_rep_mat))) y_rep_mat <- matrix(y_rep_mat, ncol = length(yrep_pattern))
    colnames(y_rep_mat) <- gsub(".*\\[|\\]", "", yrep_pattern)
    as.data.frame(y_rep_mat)
  }

  chain_dfs <- lapply(draws_list, extract_chain)
  draws_df <- dplyr::bind_rows(chain_dfs, .id = "chain") %>%
    dplyr::mutate(iter = stats::ave(rep(1, dplyr::n()), chain, FUN = seq_along))

  draws_long <- draws_df %>%
    tidyr::pivot_longer(-c(.data$chain, .data$iter), names_to = "time_index", values_to = "y_rep") %>%
    dplyr::mutate(time_index = as.numeric(.data$time_index))

  alpha <- (1 - ci) / 2
  ribbon_df <- draws_long %>%
    dplyr::group_by(.data$time_index) %>%
    dplyr::summarise(
      median = stats::median(.data$y_rep),
      lower = stats::quantile(.data$y_rep, alpha),
      upper = stats::quantile(.data$y_rep, 1 - alpha),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$time_index)

  if (!is.null(data) && "time" %in% names(data)) {
    time_mapping <- data.frame(time_index = seq_len(nrow(data)), time = data$time)
    ribbon_df <- ribbon_df %>%
      dplyr::left_join(time_mapping, by = "time_index") %>%
      dplyr::filter(!is.na(.data$time))
  } else if (!is.null(fit$time)) {
    time_mapping <- data.frame(time_index = seq_along(fit$time), time = fit$time)
    ribbon_df <- ribbon_df %>%
      dplyr::left_join(time_mapping, by = "time_index") %>%
      dplyr::filter(!is.na(.data$time))
  } else {
    ribbon_df <- ribbon_df %>%
      dplyr::mutate(time = .data$time_index)
  }

  p <- ggplot2::ggplot(ribbon_df, ggplot2::aes(x = .data$time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), fill = "gray80", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = .data$median), color = "black") +
    ggplot2::labs(x = "Time", y = "Predicted Count") +
    ggplot2::theme_bw()

  if (!is.null(data)) {
    stopifnot(all(c("time", "count") %in% names(data)))
    p <- p + ggplot2::geom_point(data = data, ggplot2::aes(x = .data$time, y = .data$count), color = "black")
  }

  if (!is.null(shadow_breakpoints)) {
    p = add_breakpoint_shadows(p, shadow_breakpoints, colors = shadow_colors)
  }

  p
}

#' Plot parameter posterior distributions
#'
#' Plots posterior densities for specified parameters from MCMC draws.
#'
#' @param fit_draws A list of chains, where each chain is a named list of parameter draws.
#' @param params Character vector of parameters (or regex patterns) to extract.
#' @param faceted Logical; whether to facet the plot by parameter (default `TRUE`).
#'
#' @return A `ggplot` object showing posterior densities for the selected parameters.
#' @examples
#' \dontrun{
#'   plot_parameter_posteriors(fit$draws, params = c("rho[1]", "K"))
#' }
#'
#' @export
plot_parameter_posteriors <- function(fit_draws, params, faceted = TRUE, colors = NULL) {

  chains <- fit_draws
  if (!is.list(chains)) stop("fit_draws must be a list of chains.")

  extract_params <- function(chain) {
    param_names <- names(chain)
    matched <- unique(unlist(lapply(params, function(p) if (p %in% param_names) p else grep(p, param_names, value = TRUE))))
    if (length(matched) == 0) stop(sprintf("No parameters matching '%s' found.", paste(params, collapse = ", ")))
    param_mat <- sapply(matched, function(p) unlist(chain[[p]]))
    if (is.null(dim(param_mat))) param_mat <- matrix(param_mat, ncol = length(matched))
    colnames(param_mat) <- matched
    as.data.frame(param_mat)
  }

  chain_dfs <- lapply(chains, extract_params)
  draws_df <- dplyr::bind_rows(chain_dfs, .id = "chain") %>%
    dplyr::mutate(iter = stats::ave(rep(1, n()), chain, FUN = seq_along))

  draws_long <- draws_df %>% tidyr::pivot_longer(-c(chain, iter), names_to = "parameter", values_to = "value")

  if (is.null(colors)) {
    colors = get_group_colors()
  }

  p <- ggplot2::ggplot(draws_long, ggplot2::aes(x = value, fill = parameter)) +
    ggplot2::geom_density(alpha = 0.7) +
    ggplot2::labs(x = "Value", y = "Density", fill = "Parameter") +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values = colors)

  if (faceted) {
    p <- p + ggplot2::facet_wrap(~parameter, scales = "free", ncol = 2) +
      ggplot2::theme(legend.position = "none")
  }
  p
}

#' Plot MCMC trace plots for parameters
#'
#' Generates trace plots for parameters in an MCMC fit.
#'
#' @param fit A fitted model object containing posterior draws in `fit$draws`.
#'
#' @return A trace plot generated by `bayesplot::mcmc_trace`.
#' @examples
#' \dontrun{
#'   plot_trace_parameters(fit)
#' }
#'
#' @export
plot_trace_parameters <- function(fit) {
  draws_df <- posterior::as_draws_df(fit$draws, parameters)
  bayesplot::mcmc_trace(draws_df)
}
