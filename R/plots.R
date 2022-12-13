
COLORS <- c("darkgreen", "dodgerblue", "darkorange", "darkorchid3")

#' Plot the evolution of the population over time
#'
#' @param d Data frame. Must have at least three columns with the following names:
#' "time", "count" or "pop.size", and "sample".
#' @param estimate_rates Boolean. If TRUE, curves computed by a rough estimation of the
#' growth rates for each time interval will be plotted.
#' @returns A plot. Represents the evolution of the population over time.
#' @import ggplot2
#' @importFrom dplyr tibble
#' @importFrom glue glue
#' @export
evolution_plot = function(d, estimate_rates = FALSE) {
  check_input(d)
  d = clean_input(d)

  sample = d$sample[1]

  if (estimate_rates) {
    delta_t = d$time[2:nrow(d)] - d$time[1:nrow(d)-1]
    n_ratio = d$count[2:nrow(d)] / d$count[1:nrow(d)-1]
    lambda = 1 / delta_t * log(n_ratio)

    pop.plot <- ggplot()

    for (i in 2:nrow(d)) {
      dt = delta_t[i-1]
      new_xs = seq(0, dt, by = 0.01)
      new_ys = d$count[i-1] * exp(new_xs * lambda[i-1])
      new_xs = new_xs + d$time[i-1]
      col = c(col, rep(COLORS[i-1], length(new_xs)))

      D = tibble(time=new_xs, count = new_ys)
      pop.plot <- pop.plot + geom_line(data=D, aes(x=time, y=count), show.legend = F, col=COLORS[i-1])
    }

    pop.plot <- pop.plot +
      geom_point(data=d, aes(x=time, y=count), col="black") +
      ggtitle(glue("{sample}")) +
      theme_bw()
    return(pop.plot)
  }

  ggplot() +
    geom_line(data=d, aes(x=time, y=count), col=COLORS[1:nrow(d)]) +
    geom_point(data=d, aes(x=time, y=count), col="black") +
    ggtitle(glue("{sample}")) +
    theme_bw() -> pop.plot
  pop.plot
}

#' Plot the posterior distribution for the birth rates.
#'
#' @param fit_model Is the result obtained by fitting a stan model over the data.
#' @returns A plot. Represents the posterior distribution for the birth rates.
#' @import ggplot2
#' @export
birth_posterior_plot = function(fit_model) {
  # plot density of growth
  posterior <- as.data.frame(fit_model)

  posterior.gathered <- tidyr::gather(posterior)
  lambda.posterior <- posterior.gathered[posterior.gathered$key %in% unique(posterior.gathered$key)[grep("^lambda", unique(posterior.gathered$key))],]
  max_x <- max(as.numeric(unlist(posterior.gathered[!(posterior.gathered$key %in% c("lp__")),2])))

  n_params = length(unique(lambda.posterior$key))
  if (n_params == 1) color_values = "gray" else color_values = COLORS[1:n_params]

  lambda.plot <- ggplot(lambda.posterior, aes(x=value, fill=key)) +
    geom_density(alpha=.5) +
    theme_bw() +
    scale_fill_manual(
      values = color_values,
      name = "Birth rates") +
    ggtitle("Posterior density for birth rates") +
    theme(plot.title = element_text(size=16),
          axis.title = element_text(size=14)) +
    xlim(0, max_x)
  lambda.plot
}

#' Plot the posterior distribution for the birth rates.
#'
#' @param fit_model Is the result obtained by fitting a stan model over the data.
#' @returns A plot. Represents the posterior distribution for the death rates.
#' @import ggplot2
#' @export
death_posterior_plot = function(fit_model) {
  # plot posteriors for deltas
  posterior <- as.data.frame(fit_model)

  posterior.gathered <- tidyr::gather(posterior)
  mu.posterior <- posterior.gathered[posterior.gathered$key %in% unique(posterior.gathered$key)[grep("^mu", unique(posterior.gathered$key))],]
  max_x <- max(as.numeric(unlist(posterior.gathered[!(posterior.gathered$key %in% c("lp__")),2])))

  n_params = length(unique(mu.posterior$key))
  if (n_params == 1) color_values = "gray" else color_values = COLORS[1:n_params]

  mu.plot <- ggplot(mu.posterior, aes(x=value, fill=key)) +
    geom_density(alpha=.5) +
    theme_bw() +
    scale_fill_manual(
      values = color_values,
      name = "Death rates") +
    ggtitle("Posterior density for death rates") +
    theme(plot.title = element_text(size=16),
          axis.title = element_text(size=14)) +
    xlim(0, max_x)
  mu.plot
}
