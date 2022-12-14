
COLORS <- c("darkgreen", "dodgerblue", "darkorange", "darkorchid3")
base_color <- "IndianRed3"

#' Plot the evolution of the population over time
#'
#' @param x A biPOD object of class `bipod`.
#' @param estimate_rates Boolean. If TRUE, curves computed by a rough estimation of the
#' growth rates for each time interval will be plotted.
#' @returns A plot. Represents the evolution of the population over time.
#' @import ggplot2
#' @export
evolution_plot <- function(x, estimate_rates = FALSE) {
  stopifnot(inherits(x, 'bipod'))

  sample <- x$sample
  times <- x$counts$time
  counts <- x$counts$count
  l <- length(times)

  if (estimate_rates) {
    delta_t <- times[2:l] - times[1:l - 1]
    n_ratio <- counts[2:l] / counts[1:l - 1]
    lambda <- 1 / delta_t * log(n_ratio)

    pop.plot <- ggplot()

    for (i in 2:length(times)) {
      dt <- delta_t[i - 1]
      new_xs <- seq(0, dt, by = 0.01)
      new_ys <- counts[i - 1] * exp(new_xs * lambda[i - 1])
      new_xs <- new_xs + times[i - 1]
      col <- c(col, rep(COLORS[i - 1], length(new_xs)))

      D <- tibble(time = new_xs, count = new_ys)
      pop.plot <- pop.plot + geom_line(data = D, aes(x = time, y = count), show.legend = F, col = COLORS[i - 1])
    }

    pop.plot <- pop.plot +
      geom_point(data = x$counts, aes(x = time, y = count), col = "black") +
      ggtitle(sample) +
      theme_bw()
    return(pop.plot)
  }

  ggplot() +
    geom_line(data = x$counts, aes(x = time, y = count), col = COLORS[1:nrow(x$counts)]) +
    geom_point(data = x$counts, aes(x = time, y = count), col = "black") +
    ggtitle(sample) +
    theme_bw() -> pop.plot
  pop.plot
}

#' Plot the posterior distribution for a specific group of parameters.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param param_name A string indicating the desired parameters.
#' @returns A plot. Represents the posterior distribution for the desired parameters.
#' @import ggplot2
#' @importFrom tidyr gather
#' @export
posterior_plot <- function(x, param_name) {
  stopifnot(inherits(x, "bipod"))
  stopifnot(param_name %in% c("lambda", "mu", "birth", "death"))

  if (param_name %in% c("lambda", "birth")) {
    rex = "^lambda"
    legend_name = "Birth rates"
    title = "Posterior density for birth rates"
  } else {
    rex = "^mu"
    legend_name = "Death rates"
    title = "Posterior density for death rates"
  }

  # plot density of growth
  fit_model <- x$fit
  posterior <- as.data.frame(fit_model)

  posterior.gathered <- tidyr::gather(posterior)
  pos <- posterior.gathered[posterior.gathered$key %in% unique(posterior.gathered$key)[grep(rex, unique(posterior.gathered$key))], ]
  max_x <- max(pos$value)

  n_params <- length(unique(pos$key))
  if (n_params == 1) color_values <- base_color else color_values <- COLORS[1:n_params]

  p <- ggplot(pos, aes(x = value, fill = key)) +
    geom_density(alpha = .5) +
    theme_bw() +
    scale_fill_manual(
      values = color_values,
      name = legend_name
    ) +
    ggtitle(title) +
    theme(
      plot.title = element_text(size = 16),
      axis.title = element_text(size = 14)
    ) +
    xlim(0, max_x)
  p
}


#' Plot the posterior predictive checks for the population counts.
#'
#' @param x A biPOD object of class `bipod`. Must contains 'fit'
#' @param prob A value between `0` and `1` indicating the desired probability
#'   mass to include in the intervals of the generated samples.
#' @returns A plot. Represents the posterior predictive checks for the population counts.
#' @import ggplot2
#' @export
ppc_plot <- function(x, prob) {
  stopifnot(inherits(x, "bipod"))
  stopifnot((prob <= 1 && prob >= 0))

  y_data <- .get_ppc_data(x, prob)

  p <- ggplot(y_data) +
    geom_bar(
      aes(x=t, y=y, fill=y_col),
      stat = "identity",
      alpha = .5
    ) +
    geom_pointrange(
      aes(x=t, ymin=l, ymax=h, y=m, colour=y_rep_col),
      size = 1,
      fatten = 1.5
    ) +
    scale_color_manual(values = base_color) +
    scale_fill_manual(values = base_color) +
    labs(fill = "", colour="") +
    ggtitle("PPC") +
    theme(
      plot.title = element_text(size = 16),
      axis.title = element_text(size = 14)
    ) +
    scale_x_continuous(breaks = pretty) +
    theme_bw()
  p
}

.get_ppc_data = function(x, prob) {
  # extract fit and data
  fit <- x$fit

  t <- x$counts$time
  y <- x$counts$count
  n = length(y)
  t <- t[2:n]
  y <- y[2:n]

  y_generated <- as.data.frame(rstan::extract(fit, pars="N_rep")$N_rep)

  # prepare data
  y_data <- data.frame(t = t, y=y)

  alpha <- (1 - prob) / 2
  probs <- sort(c(alpha, 0.5, 1 - alpha))

  # Prepare for final summary
  lo  <- function(x) quantile(x, probs[1])
  mid <- function(x) quantile(x, probs[2])
  hi  <- function(x) quantile(x, probs[3])
  summary_funs <- list(l = lo, m = mid, h = hi)
  summary_names <- c("l", "m", "h")

  y_generated %>%
    dplyr::summarise_all(summary_funs) -> y_gen_summary

  for (sn in summary_names) {
    y_gen_summary %>%
      dplyr::select(contains(sn)) -> v
    y_data[sn] = as.numeric(v[1,])
  }
  y_data <- y_data %>%
    dplyr::mutate(y_rep_col = "y_rep", y_col = "y")
  y_data
}
