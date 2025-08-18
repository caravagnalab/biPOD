#' Plot posterior distributions of inferred breakpoints
#'
#' Plots posterior densities of the inferred breakpoint times (`t_array[j]`)
#' from a fitted breakpoint model.
#'
#' @param bp_fit A list containing the fitted breakpoint model, with element
#' `final_fit$draws` as a list of chains (each chain is a named list of parameter draws).
#' @param data Optional data frame with column `time` to restrict the x-axis to observed times.
#'
#' @return A `ggplot` object showing the posterior densities for each breakpoint parameter.
#'
#' @export
plot_breakpoint_posterior <- function(bp_fit, data = NULL, colors = NULL) {

  chains <- bp_fit$final_fit$draws
  if (!is.list(chains)) stop("bp_fit$final_fit$draws must be a list of chains.")

  extract_breakpoints <- function(chain) {
    param_names <- names(chain)
    t_array_pattern <- grep("t_array\\[", param_names, value = TRUE)
    if (length(t_array_pattern) == 0) stop("No t_array[...] parameters found in the draws.")
    t_array_mat <- sapply(t_array_pattern, function(p) unlist(chain[[p]]))
    if (is.null(dim(t_array_mat))) t_array_mat <- matrix(t_array_mat, ncol = length(t_array_pattern))
    colnames(t_array_mat) <- gsub(".*\\[|\\]", "", t_array_pattern)
    as.data.frame(t_array_mat)
  }

  chain_dfs <- lapply(chains, extract_breakpoints)
  draws_df <- dplyr::bind_rows(chain_dfs, .id = "chain") %>%
    dplyr::mutate(iter = stats::ave(rep(1, dplyr::n()), chain, FUN = seq_along))

  draws_long <- draws_df %>%
    tidyr::pivot_longer(-c(.data$chain, .data$iter), names_to = "breakpoint_index", values_to = "time_value") %>%
    dplyr::mutate(
      breakpoint_index = as.numeric(.data$breakpoint_index),
      breakpoint_label = paste0("t_array[", .data$breakpoint_index, "]")
    )

  if (is.null(colors)) {
    colors = get_group_colors()
  }

  p <- ggplot2::ggplot(draws_long, ggplot2::aes(x = time_value, fill = breakpoint_label)) +
    ggplot2::geom_density(alpha = 0.7) +
    ggplot2::labs(x = "Time", y = "Density", fill = "Breakpoint") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_fill_manual(values = colors)

  if (!is.null(data) && "time" %in% names(data)) {
    p <- p + ggplot2::scale_x_continuous(limits = c(min(data$time), max(data$time)))
  }

  p
}

#' Plot model selection results for breakpoint detection
#'
#' Visualizes the model comparison scores (BIC or LOOIC) across candidate segmentations,
#' highlighting the best number of segments.
#'
#' @param bp_fit A list returned by \code{\link{fit_breakpoints}} containing
#'   `evaluation_table` with columns `num_segments`, `score`, and `criterion`,
#'   as well as `final_breakpoints` to identify the best segmentation.
#'
#' @return A `ggplot` object showing the model selection curve with the best
#' segmentation highlighted.
#'
#' @export
plot_breakpoint_model_selection <- function(bp_fit) {

  best_num_segments <- length(bp_fit$final_breakpoints) + 1
  df <- bp_fit$evaluation_table %>%
    dplyr::mutate(is_best = .data$num_segments == best_num_segments)

  criterion <- unique(df$criterion)
  y_label <- ifelse(criterion == "bic", "BIC", "LOOIC")

  ggplot2::ggplot(df, ggplot2::aes(x = .data$num_segments, y = .data$score)) +
    ggplot2::geom_line(color = "grey60") +
    ggplot2::geom_point(size = 3, color = "black") +
    ggplot2::geom_point(
      data = dplyr::filter(df, .data$is_best),
      ggplot2::aes(x = .data$num_segments, y = .data$score),
      size = 5, shape = 1, color = "red", stroke = 1.2
    ) +
    ggplot2::labs(
      title = "Model Selection for Breakpoint Detection",
      x = "Number of Segments",
      y = y_label
    ) +
    ggplot2::theme_bw()
}
