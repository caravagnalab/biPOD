#' Fit the joint clone-selection model
#'
#' Joint piecewise-exponential model of the read counts supporting each clone's
#' markers (for example clone-specific structural variants in cell-free DNA).
#' All clones share one breakpoint `t_b` (response to regrowth); each clone `k`
#' has its own initial size, response-phase rate `rho1[k]` and post-breakpoint
#' rate `rho2[k]`. Counts are negative binomial with a per-sample depth offset.
#' Unlike the two-population recovery models, the initial size of each clone is
#' free, so clones that are small or undetected at baseline are allowed.
#'
#' Selection is summarised by `delta_log_freq`, the change in each clone's log
#' relative frequency over the window, and `p_selected`, its posterior
#' probability of being positive.
#'
#' Chains start with `t_b` at its prior centre. Setting a small `t_b_prior_sd`
#' effectively fixes the breakpoint (for example at a regrowth time estimated
#' from a denser marker). When a clone grows in a held phase, `rho2` averages
#' the hold and the regrowth, so it underestimates the regrowth rate.
#'
#' @param data Long data frame with columns `time`, `clone` and `count` (reads
#'   supporting the clone's markers), and optionally `depth` (total reads over
#'   the same markers). Without `depth`, all samples get the same depth.
#' @param t_b_prior Prior centre of the shared breakpoint. Defaults to the time
#'   of the lowest total depth-normalised count.
#' @param t_b_prior_sd Prior sd of the breakpoint. Defaults to 10% of the time span.
#' @param rho_prior_sd Prior sd of the clone rates.
#' @param chains Number of MCMC chains.
#' @param iter Number of warmup and of sampling iterations per chain.
#' @param seed Random seed.
#' @param cores Number of CPU cores.
#' @param adapt_delta Target acceptance rate for NUTS.
#' @param max_rhat Rhat threshold used for the `converged` flag.
#'
#' @return A list containing:
#'   \item{clones}{Tibble with one row per clone: posterior median and 90%
#'     interval of `rho1`, `rho2` and `delta_log_freq`, `p_selected`, and the
#'     number of non-finite `delta_log_freq` draws (dropped from the summary).}
#'   \item{t_b}{Posterior median and 90% interval of the breakpoint.}
#'   \item{fit}{Parsed Stan fit.}
#'   \item{max_rhat}{Maximum Rhat over the sampled parameters.}
#'   \item{n_divergent}{Number of divergent transitions.}
#'   \item{converged}{Whether `max_rhat` is below the threshold.}
#'
#' @examples
#' \dontrun{
#'   d <- data.frame(time = rep(c(0, 30, 60, 200, 400), 2),
#'                   clone = rep(c("A", "B"), each = 5),
#'                   count = c(500, 40, 3, 0, 0, 20, 2, 1, 15, 300),
#'                   depth = 1000)
#'   fit_clone_selection(d)
#' }
#' @export
fit_clone_selection <- function(data,
                                t_b_prior = NULL,
                                t_b_prior_sd = NULL,
                                rho_prior_sd = 0.1,
                                chains = 4,
                                iter = 2000,
                                seed = 123,
                                cores = 4,
                                adapt_delta = 0.95,
                                max_rhat = 1.05) {
  stopifnot(all(c("time", "clone", "count") %in% colnames(data)))
  if (!("depth" %in% colnames(data))) data$depth <- 1
  if (any(data$count < 0)) stop("'count' must be non-negative")

  # clones without a single supporting read carry no information
  keep <- tapply(data$count, data$clone, function(x) any(x > 0))
  data <- data[data$clone %in% names(keep)[keep], ]
  clones <- sort(unique(as.character(data$clone)))
  if (length(clones) < 2) stop("At least two clones with supporting reads are needed")
  times <- sort(unique(data$time))
  if (length(times) < 3) stop("At least three time points are needed")

  y <- matrix(0L, length(clones), length(times))
  off <- matrix(0, length(clones), length(times))
  for (k in seq_along(clones)) {
    dk <- data[data$clone == clones[k], ]
    idx <- match(dk$time, times)
    if (anyNA(idx) || length(idx) != length(times)) {
      stop(sprintf("Clone '%s' must have exactly one row per time point", clones[k]))
    }
    y[k, idx] <- as.integer(round(dk$count))
    off[k, idx] <- log(dk$depth / stats::median(dk$depth))
  }

  if (is.null(t_b_prior)) {
    t_b_prior <- times[which.min(colSums(y / exp(off)))]
  }
  if (is.null(t_b_prior_sd)) t_b_prior_sd <- 0.1 * (max(times) - min(times))

  stan_data <- list(
    K = length(clones), S = length(times), T = times, y = y, depth_offset = off,
    a_prior = log(y[, 1] + 0.5) - off[, 1],
    t_b_prior = t_b_prior, t_b_prior_sd = t_b_prior_sd, rho_prior_sd = rho_prior_sd
  )

  # start every chain at the breakpoint prior (needed when the prior is narrow)
  t_b0 <- min(max(t_b_prior, min(times) + 1e-6), max(times) - 1e-6)
  inits <- lapply(seq_len(chains), function(i) list(t_b = t_b0))
  mod <- get_model("clone_selection", "negbinomial")
  out_dir <- tempfile("biPOD_clones_")
  dir.create(out_dir)
  fit <- suppressMessages(suppressWarnings(mod$sample(
    data = stan_data, chains = chains, iter_warmup = iter, iter_sampling = iter,
    seed = seed, parallel_chains = cores, refresh = 0, init = inits,
    adapt_delta = adapt_delta, output_dir = out_dir, show_messages = FALSE
  )))

  draws <- posterior::as_draws_df(fit$draws(c("rho1", "rho2", "delta_log_freq", "t_b")))
  q <- function(x) stats::quantile(x, c(0.5, 0.05, 0.95), names = FALSE)
  clone_tbl <- dplyr::bind_rows(lapply(seq_along(clones), function(k) {
    col <- function(v) {
      x <- draws[[sprintf("%s[%d]", v, k)]]
      x[is.finite(x)]
    }
    r1 <- q(col("rho1")); r2 <- q(col("rho2")); dl <- col("delta_log_freq"); dq <- q(dl)
    dplyr::tibble(
      clone = clones[k],
      rho1 = r1[1], rho1_q5 = r1[2], rho1_q95 = r1[3],
      rho2 = r2[1], rho2_q5 = r2[2], rho2_q95 = r2[3],
      delta_log_freq = dq[1], delta_q5 = dq[2], delta_q95 = dq[3],
      p_selected = mean(dl > 0),
      n_nonfinite_delta = sum(!is.finite(draws[[sprintf("delta_log_freq[%d]", k)]]))
    )
  }))

  rhat <- max(fit$summary(c("a", "rho1", "rho2", "t_b", "phi"))$rhat, na.rm = TRUE)
  n_div <- sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
  if (rhat >= max_rhat) {
    cli::cli_alert_warning("Clone-selection model did not converge (max Rhat {round(rhat, 2)}).")
  }

  list(clones = clone_tbl, t_b = stats::setNames(q(draws$t_b), c("median", "q5", "q95")),
       fit = parse_stan_fit(fit), max_rhat = rhat, n_divergent = n_div, converged = rhat < max_rhat)
}
