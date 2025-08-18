
#' Plot posterior predictive ribbon for U-shape dynamic
#'
#' Creates a ribbon plot showing the median and credible intervals of the posterior.
#' It plots both total population as well as resistant and sensitive (if detected)
#'
#' @param fit .
#' @param data Optional data frame with columns `time` and `count` to overlay observed data.
#' @param ci Credible interval width (default `0.9`).
#'
#' @return A `ggplot` object showing the posterior predictive ribbon and observed data.
#'
#' @export
plot_u_ribbon <- function(fit, data = NULL, ci = 0.9) {

  draws_list <- fit$final_fit$draws
  if (!is.list(draws_list)) stop("fit$final_fit$draws must be a list of chains.")

  # Extract y_rep / ns / nr draws from one chain into a data.frame of iterations x parameters
  extract_chain <- function(chain) {
    param_names <- names(chain)

    # Accept y_rep[i] or yrep[i], plus ns[i] and nr[i]
    keep <- grep("^(y[_]?rep|ns|nr)\\[\\d+\\]$", param_names, value = TRUE)
    if (length(keep) == 0) {
      stop("No parameters found among y_rep[i]/yrep[i], ns[i], nr[i].")
    }

    mat <- sapply(keep, function(p) unlist(chain[[p]]))
    if (is.null(dim(mat))) mat <- matrix(mat, ncol = length(keep))
    colnames(mat) <- keep
    as.data.frame(mat)
  }

  # Bind chains and tag iteration per chain
  chain_dfs <- lapply(draws_list, extract_chain)
  draws_df <- dplyr::bind_rows(chain_dfs, .id = "chain") %>%
    dplyr::mutate(iter = stats::ave(rep(1, dplyr::n()), chain, FUN = seq_along))

  # Long format: parse series (y_rep/ns/nr) and time_index
  draws_long <- draws_df %>%
    tidyr::pivot_longer(-c(.data$chain, .data$iter), names_to = "param", values_to = "value") %>%
    dplyr::mutate(
      series = sub("\\[.*$", "", .data$param),          # y_rep, yrep, ns, nr
      series = ifelse(series %in% c("yrep", "y_rep"), "y_rep", series),  # normalize
      time_index = as.numeric(gsub(".*\\[|\\]", "", .data$param))
    ) %>%
    dplyr::select(-.data$param)

  # Summarise ribbons per series and time index
  alpha <- (1 - ci) / 2
  ribbon_df <- draws_long %>%
    dplyr::group_by(.data$series, .data$time_index) %>%
    dplyr::summarise(
      median = stats::median(.data$value),
      lower  = stats::quantile(.data$value, alpha),
      upper  = stats::quantile(.data$value, 1 - alpha),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$series, .data$time_index)

  # Map time_index -> time (from data$time or fit$time, else use index)
  if (!is.null(data) && "time" %in% names(data)) {
    time_mapping <- data.frame(time_index = seq_len(nrow(data)), time = data$time)
  } else if (!is.null(fit$time)) {
    time_mapping <- data.frame(time_index = seq_along(fit$time), time = fit$time)
  } else {
    time_mapping <- data.frame(time_index = sort(unique(ribbon_df$time_index)),
                               time = sort(unique(ribbon_df$time_index)))
  }

  ribbon_df <- ribbon_df %>%
    dplyr::left_join(time_mapping, by = "time_index")  %>%
    dplyr::filter(!is.na(.data$time))

  ribbon_df$series[ribbon_df$series == "nr"] = "N resistant"
  ribbon_df$series[ribbon_df$series == "ns"] = "N sensitive"
  ribbon_df$series[ribbon_df$series == "y_rep"] = "N total"

  # Plot: ribbons + medians for each series
  p <- ggplot2::ggplot(ribbon_df, ggplot2::aes(x = .data$time)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data$series),
      alpha = 0.25
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$median, color = .data$series)) +
    ggplot2::labs(x = "Time", y = "Predicted value", fill = "Series", color = "Series") +
    ggplot2::theme_bw() +
    scale_color_manual(values = c("N total" = "black", "N sensitive" = "mediumpurple", "N resistant" = "goldenrod")) +
    scale_fill_manual(values = c("N total" = "black", "N sensitive" = "mediumpurple", "N resistant" = "goldenrod"))

  # Optional observed data overlay
  if (!is.null(data)) {
    stopifnot(all(c("time", "count") %in% names(data)))
    p <- p + ggplot2::geom_point(
      data = data,
      ggplot2::aes(x = .data$time, y = .data$count),
      inherit.aes = FALSE,
      color = "black"
    )
  }

  p
}
