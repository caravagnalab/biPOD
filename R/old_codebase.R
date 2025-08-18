#
# rm(list = ls())
# library(cmdstanr)
# library(loo)
# library(posterior)
# library(ggplot2)
# library(tidyverse)
# library(bayesplot)
#
# fit_growth_model_breakpoints <- function(data,
#                                          breakpoints,
#                                          stan_dir,
#                                          chains = 4,
#                                          iter = 4000,
#                                          seed = 123,
#                                          cores = 4,
#                                          t_prior_sd = 1) {
#   # Ensure data is ordered by time
#   stopifnot(all(c("time", "count") %in% colnames(data)))
#   data <- data[order(data$time), ]
#
#   # Prepare Stan data
#   S <- nrow(data)
#   G <- length(breakpoints) + 1
#   T_obs <- data$time
#   N_obs <- data$count
#
#   stan_data <- list(
#     S = S,
#     G = G,
#     N = N_obs,
#     T = T_obs,
#     t_prior = as.vector(breakpoints),
#     t_prior_sd = t_prior_sd
#   )
#
#   # Compile model
#   model_file <- file.path(stan_dir, "exponential_with_init_bp.stan")
#   mod <- cmdstan_model(model_file)
#
#   # Fit model
#   fit <- mod$sample(
#     data = stan_data,
#     chains = chains,
#     iter_warmup = iter / 2,
#     iter_sampling = iter / 2,
#     seed = seed,
#     parallel_chains = cores,
#     refresh = 500
#   )
#
#   return(fit)
# }
#
# search_best_breakpoints <- function(data,
#                                     max_segments = 4,
#                                     min_segment_size = 3,
#                                     comparison = "bic",
#                                     ...) {
#
#   require(dplyr)
#   require(digest)
#
#   # Cache to avoid re-fitting the same configuration
#   cache <- new.env()
#
#   evaluate_model <- function(breaks) {
#     key <- digest::digest(sort(breaks))
#     if (exists(key, envir = cache)) {
#       return(get(key, envir = cache))
#     } else {
#       res <- fit_growth_models(data = data, breakpoints = breaks, comparison = comparison, ...)
#       assign(key, res, envir = cache)
#       return(res)
#     }
#   }
#
#   # Container for best model per segment count
#   best_per_k <- list()
#   candidate_times <- unique(data$time)
#   candidate_times <- candidate_times[candidate_times > min(data$time) & candidate_times < max(data$time)]
#
#   # Try each number of segments from 1 to max_segments
#   for (k in 1:max_segments) {
#     message(sprintf("Evaluating %d segment(s)", k))
#     if (k == 1) {
#       # No breakpoints
#       result <- evaluate_model(numeric(0))
#       score <- if (comparison == "bic") result$model_table$BIC else result$model_table[1, "looic"]
#       best_per_k[[1]] <- list(result = result, breakpoints = numeric(0), score = score)
#       next
#     }
#
#     # Try all combinations of k-1 breakpoints
#     combn_breaks <- combn(candidate_times, k - 1, simplify = FALSE)
#
#     best_result <- NULL
#     best_score <- Inf
#     best_breaks <- NULL
#
#     for (breaks in combn_breaks) {
#       breaks <- sort(breaks)
#
#       # Check that all segments have enough data
#       segments <- cut(data$time, breaks = c(-Inf, breaks, Inf))
#       seg_sizes <- table(segments)
#       if (any(seg_sizes < min_segment_size)) next
#
#       result <- evaluate_model(breaks)
#       score <- if (comparison == "bic") result$model_table$BIC else result$model_table[1, "looic"]
#
#       if (score < best_score) {
#         best_result <- result
#         best_score <- score
#         best_breaks <- breaks
#       }
#     }
#
#     best_per_k[[k]] <- list(result = best_result, breakpoints = best_breaks, score = best_score)
#   }
#
#   # Find best number of segments
#   scores <- sapply(best_per_k, function(x) x$score)
#   best_k <- which.min(scores)
#
#   list(
#     best_k = best_k,
#     best_model = best_per_k[[best_k]],
#     all_models = best_per_k,
#     comparison = comparison
#   )
# }
#
# fit_growth_models <- function(data, breakpoints, with_initiation = TRUE,
#                               stan_dir = "stan_models", chains = 4, iter = 2000, seed = 123, cores = 4,
#                               comparison = c("loo", "bic"),
#                               models_to_fit = c("exponential", "logistic", "gompertz")) {
#   comparison <- match.arg(comparison)
#
#   stopifnot(all(c("time", "count") %in% colnames(data)))
#   data <- data[order(data$time), ]
#
#   t_array <- if (length(breakpoints) == 0) array(0, dim = 0) else as.vector(breakpoints)
#   G <- length(breakpoints) + 1
#
#   stan_data <- list(
#     S = nrow(data),
#     G = G,
#     N = data$count,
#     T = data$time,
#     t_array = t_array
#   )
#
#   model_files <- paste0(models_to_fit, if (with_initiation) "_with_init.stan" else "_no_init.stan")
#
#   fits <- list()
#   comparisons <- list()
#   info_criteria <- numeric(length(models_to_fit))
#   names(info_criteria) <- models_to_fit
#
#   for (i in seq_along(models_to_fit)) {
#     model_name <- models_to_fit[i]
#     model_path <- file.path(stan_dir, model_files[i])
#     model <- cmdstan_model(model_path)
#
#     message(sprintf("Fitting model: %s", model_name))
#     fit <- suppressMessages(
#       suppressWarnings(
#         model$sample(
#           data = stan_data,
#           chains = chains,
#           iter_warmup = iter,
#           iter_sampling = iter,
#           seed = seed,
#           parallel_chains = cores,
#           refresh = 0
#         )
#       )
#     )
#
#     fits[[model_name]] <- fit
#     draws <- fit$draws(format = "draws_matrix")
#     log_lik <- draws[, grep("^log_lik", colnames(draws)), drop = FALSE]
#
#     if (comparison == "loo") {
#       comparisons[[model_name]] <- loo::loo(log_lik)
#       info_criteria[i] <- comparisons[[model_name]]$estimates["looic", "Estimate"]
#     } else if (comparison == "bic") {
#       log_mean_lik <- apply(log_lik, 2, mean)
#       log_lik_sum <- sum(log_mean_lik)
#
#       param_cols <- grep("^rho\\[|^n0$|^t0$|^K$", colnames(draws), value = TRUE)
#       n_params <- length(param_cols)
#
#       N <- length(data$count)
#       info_criteria[i] <- -2 * log_lik_sum + log(N) * n_params
#     }
#   }
#
#   comp_table <- if (comparison == "bic") {
#     data.frame(BIC = sort(info_criteria))
#   } else {
#     tbl <- loo::loo_compare(comparisons)
#     tbl <- as.data.frame(tbl)
#     tbl$model <- rownames(tbl)
#     tbl <- tbl[, c("model", setdiff(names(tbl), "model"))]
#     tbl
#   }
#
#   list(fits = fits, comparisons = comparisons, model_table = comp_table, criterion = comparison)
# }
#
# fit_growth_models_VI <- function(data, breakpoints, with_initiation = TRUE,
#                               stan_dir = "stan_models", chains = 4, iter = 2000, seed = 123, cores = 4,
#                               comparison = c("loo", "bic"),
#                               models_to_fit = c("exponential", "logistic", "gompertz"),
#                               method = c("sampling", "vi"),
#                               use_elbo = FALSE) {
#
#   comparison <- match.arg(comparison)
#   method <- match.arg(method)
#
#   stopifnot(all(c("time", "count") %in% colnames(data)))
#   data <- data[order(data$time), ]
#
#   t_array <- if (length(breakpoints) == 0) array(0, dim = 0) else as.vector(breakpoints)
#   G <- length(breakpoints) + 1
#
#   stan_data <- list(
#     S = nrow(data),
#     G = G,
#     N = data$count,
#     T = data$time,
#     t_array = t_array
#   )
#
#   model_files <- paste0(models_to_fit, if (with_initiation) "_with_init.stan" else "_no_init.stan")
#
#   fits <- list()
#   comparisons <- list()
#   info_criteria <- numeric(length(models_to_fit))
#   names(info_criteria) <- models_to_fit
#
#   for (i in seq_along(models_to_fit)) {
#     model_name <- models_to_fit[i]
#     model_path <- file.path(stan_dir, model_files[i])
#     model <- cmdstanr::cmdstan_model(model_path)
#
#     message(sprintf("Fitting model: %s (%s)", model_name, method))
#
#     if (method == "sampling") {
#       fit <- suppressMessages(
#         suppressWarnings(
#           model$sample(
#             data = stan_data,
#             chains = chains,
#             iter_warmup = iter,
#             iter_sampling = iter,
#             seed = seed,
#             parallel_chains = cores,
#             refresh = 0
#           )
#         )
#       )
#
#       draws <- fit$draws(format = "draws_matrix")
#       log_lik <- draws[, grep("^log_lik", colnames(draws)), drop = FALSE]
#
#       if (comparison == "loo") {
#         comparisons[[model_name]] <- loo::loo(log_lik)
#         info_criteria[i] <- comparisons[[model_name]]$estimates["looic", "Estimate"]
#       } else if (comparison == "bic") {
#         log_mean_lik <- apply(log_lik, 2, mean)
#         log_lik_sum <- sum(log_mean_lik)
#         param_cols <- grep("^rho\\[|^n0$|^t0$|^K$", colnames(draws), value = TRUE)
#         n_params <- length(param_cols)
#         N <- length(data$count)
#         info_criteria[i] <- -2 * log_lik_sum + log(N) * n_params
#       }
#
#     } else if (method == "vi") {
#       fit <- suppressMessages(
#         suppressWarnings(
#           model$variational(
#             data = stan_data,
#             seed = seed,
#             output_samples = iter
#           )
#         )
#       )
#
#       draws <- fit$draws(format = "draws_matrix")
#       log_lik_idx <- grep("^log_lik", colnames(draws))
#
#       if (use_elbo) {
#         elbo <- fit$metadata()$elbo
#         info_criteria[i] <- -elbo
#         comparisons[[model_name]] <- list(elbo = elbo)
#       } else if (length(log_lik_idx) > 0) {
#         log_lik <- draws[, log_lik_idx, drop = FALSE]
#
#         if (comparison == "loo") {
#           comparisons[[model_name]] <- loo::loo(log_lik)
#           info_criteria[i] <- comparisons[[model_name]]$estimates["looic", "Estimate"]
#         } else if (comparison == "bic") {
#           log_mean_lik <- apply(log_lik, 2, mean)
#           log_lik_sum <- sum(log_mean_lik)
#           param_cols <- grep("^rho\\[|^n0$|^t0$|^K$", colnames(draws), value = TRUE)
#           n_params <- length(param_cols)
#           N <- length(data$count)
#           info_criteria[i] <- -2 * log_lik_sum + log(N) * n_params
#         }
#       } else {
#         warning(sprintf("No log_lik found in draws for model %s — cannot compute %s", model_name, comparison))
#         info_criteria[i] <- NA
#         comparisons[[model_name]] <- NULL
#       }
#     }
#
#     fits[[model_name]] <- fit
#   }
#
#   comp_table <- if (use_elbo) {
#     data.frame(ELBO = sort(info_criteria, na.last = TRUE))
#   } else if (comparison == "bic") {
#     data.frame(BIC = sort(info_criteria, na.last = TRUE))
#   } else if (comparison == "loo") {
#     valid_comparisons <- comparisons[!vapply(comparisons, is.null, logical(1))]
#     tbl <- loo::loo_compare(valid_comparisons)
#     tbl <- as.data.frame(tbl)
#     tbl$model <- rownames(tbl)
#     tbl <- tbl[, c("model", setdiff(names(tbl), "model"))]
#     tbl
#   } else {
#     NULL
#   }
#
#   list(
#     fits = fits,
#     comparisons = comparisons,
#     model_table = comp_table,
#     criterion = if (use_elbo) "elbo" else comparison,
#     method = method
#   )
# }
#
# binary_segment_fit <- function(data,
#                                max_segments = 4,
#                                min_segment_size = 3,
#                                beam_width = 3,
#                                comparison = "bic",
#                                enforce_rho_separation = TRUE,
#                                alpha_rho = 0.05,
#                                ...) {
#   require(dplyr)
#   require(posterior)
#
#   rhos_distinct <- function(draws, alpha = 0.05) {
#     rho_idx <- grep("^rho\\[", colnames(draws))
#     rho_mat <- draws[, rho_idx, drop = FALSE]
#     G <- ncol(rho_mat)
#     if (G <= 1) return(TRUE)
#     ci_bounds <- apply(rho_mat, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
#     for (g in seq_len(G - 1)) {
#       if (!(ci_bounds[2, g] < ci_bounds[1, g + 1] || ci_bounds[2, g + 1] < ci_bounds[1, g]))
#         return(FALSE)
#     }
#     return(TRUE)
#   }
#
#   evaluate_model <- function(breaks) {
#     fit_growth_models(data = data,
#                       breakpoints = breaks,
#                       comparison = comparison,
#                       ...)
#   }
#
#   validate_breakpoints <- function(breaks) {
#     if (length(breaks) == 0) return(TRUE)
#     segments <- sort(c(min(data$time), breaks, max(data$time)))
#     for (i in seq_len(length(segments) - 1)) {
#       segment_start <- segments[i]
#       segment_end <- segments[i + 1]
#       if (sum(data$time >= segment_start & data$time < segment_end) < min_segment_size)
#         return(FALSE)
#     }
#     return(TRUE)
#   }
#
#   all_evaluations <- list()
#   all_checked_candidates <- list()
#
#   store_evaluation <- function(breaks, result) {
#     metric_col <- if (comparison == "bic") "BIC" else "looic"
#     for (growth_type in rownames(result$model_table)) {
#       all_evaluations <<- append(all_evaluations, list(data.frame(
#         num_breakpoints = length(breaks),
#         num_segments = length(breaks) + 1,
#         breakpoints = if (length(breaks) == 0) "none" else paste(round(breaks, 3), collapse = ", "),
#         growth_type = growth_type,
#         metric_value = result$model_table[growth_type, metric_col],
#         stringsAsFactors = FALSE
#       )))
#     }
#   }
#
#   # Initial fit (no breakpoints)
#   initial_result <- evaluate_model(numeric(0))
#   initial_metric <- if (comparison == "bic") min(initial_result$model_table$BIC) else initial_result$model_table[1, "looic"]
#   store_evaluation(numeric(0), initial_result)
#   candidates <- list(list(breakpoints = numeric(0), result = initial_result, metric = initial_metric))
#
#   all_checked_candidates <- append(all_checked_candidates, list(list(
#     breakpoints = numeric(0),
#     num_segments = 1,
#     valid = TRUE,
#     rejected_reason = NA,
#     result = initial_result,
#     metric = initial_metric
#   )))
#
#   for (k in 1:(max_segments - 1)) {
#     new_candidates <- list()
#
#     for (candidate in candidates) {
#       current_breaks <- candidate$breakpoints
#       segments <- sort(c(min(data$time), current_breaks, max(data$time)))
#
#       for (i in seq_len(length(segments) - 1)) {
#         segment_start <- segments[i]
#         segment_end <- segments[i + 1]
#         candidate_times <- data %>%
#           filter(time > segment_start, time < segment_end) %>%
#           pull(time) %>%
#           unique() %>%
#           sort()
#
#         if (length(candidate_times) < min_segment_size) next
#
#         for (t in candidate_times) {
#           new_breaks <- sort(unique(c(current_breaks, t)))
#           if (length(new_breaks) + 1 > max_segments) next
#
#           candidate_record <- list(
#             breakpoints = new_breaks,
#             num_segments = length(new_breaks) + 1,
#             valid = TRUE,
#             rejected_reason = NA,
#             result = NULL,
#             metric = NA
#           )
#
#           if (!validate_breakpoints(new_breaks)) {
#             candidate_record$valid <- FALSE
#             candidate_record$rejected_reason <- "min_segment_size"
#             all_checked_candidates <- append(all_checked_candidates, list(candidate_record))
#             next
#           }
#
#           result <- evaluate_model(new_breaks)
#           candidate_record$result <- result
#           store_evaluation(new_breaks, result)
#
#           if (enforce_rho_separation) {
#             model_draws <- posterior::as_draws_matrix(result$fits[[1]]$draws())
#             if (!rhos_distinct(model_draws, alpha = alpha_rho)) {
#               candidate_record$valid <- FALSE
#               candidate_record$rejected_reason <- "rho_overlap"
#               all_checked_candidates <- append(all_checked_candidates, list(candidate_record))
#               next
#             }
#           }
#
#           metric <- if (comparison == "bic") min(result$model_table$BIC) else result$model_table[1, "looic"]
#           candidate_record$metric <- metric
#
#           new_candidates <- append(new_candidates, list(list(
#             breakpoints = new_breaks,
#             result = result,
#             metric = metric
#           )))
#
#           all_checked_candidates <- append(all_checked_candidates, list(candidate_record))
#         }
#       }
#     }
#
#     if (length(new_candidates) == 0) break
#     sorted <- new_candidates[order(sapply(new_candidates, function(x) x$metric))]
#     candidates <- sorted[seq_len(min(beam_width, length(sorted)))]
#
#     invisible(lapply(candidates, function(c) {
#       message("Retained breakpoints: ", paste(round(c$breakpoints, 3), collapse = ", "))
#     }))
#   }
#
#   valid_evals <- Filter(function(x) isTRUE(x$valid) && !is.na(x$metric), all_checked_candidates)
#   if (length(valid_evals) == 0) stop("No valid candidate models were found.")
#
#   best_candidate <- valid_evals[[which.min(sapply(valid_evals, function(x) x$metric))]]
#   best_candidate$comparison_table <- do.call(rbind, all_evaluations)
#   best_candidate$comparison_metric <- comparison
#   best_candidate$all_checked_candidates <- all_checked_candidates
#
#   return(best_candidate)
# }
#
# fit_piecewise_bayes <- function(data, max_K = 3, comparison = c("bic", "loo"),
#                                 stan_file = "piecewise_linear_K.stan",
#                                 iter = 2000, chains = 4, seed = 123) {
#   comparison <- match.arg(comparison)
#   stopifnot(all(c("time", "count") %in% names(data)))
#
#   x <- data$time
#   y <- log1p(data$count)
#   N <- length(x)
#
#   fit_list <- list()
#   score_list <- numeric(max_K + 1)
#   tau_list <- list()
#
#   for (K in 0:max_K) {
#     stan_data <- list(
#       N = N,
#       x = x,
#       y = y,
#       K = K,
#       min_x = min(x),
#       max_x = max(x)
#     )
#
#     mod <- cmdstan_model(stan_file)
#     fit <- mod$sample(data = stan_data, chains = chains, iter_sampling = iter,
#                       iter_warmup = iter / 2, seed = seed, refresh = 0)
#
#     fit_list[[K + 1]] <- fit
#     draws <- fit$draws(variables = "log_lik")
#     log_lik_mat <- posterior::as_draws_matrix(draws)
#
#     if (comparison == "loo") {
#       score_list[K + 1] <- loo(log_lik_mat)$estimates["elpd_loo", "Estimate"]
#     } else if (comparison == "bic") {
#       log_lik_sum <- mean(rowSums(log_lik_mat))
#       n_params <- 1 + (K + 1) + K + 1  # alpha + beta + tau + sigma
#       score_list[K + 1] <- -2 * log_lik_sum + n_params * log(N)
#     }
#
#     # Extract median taus if any
#     if (K > 0) {
#       tau_draws <- fit$draws(variables = "tau")
#       tau_median <- apply(posterior::as_draws_matrix(tau_draws), 2, median)
#       tau_list[[K + 1]] <- sort(tau_median)
#     } else {
#       tau_list[[K + 1]] <- numeric(0)
#     }
#   }
#
#   # Select best model
#   best_idx <- if (comparison == "loo") which.max(score_list) else which.min(score_list)
#
#   return(list(
#     best_K = best_idx - 1,
#     best_score = score_list[best_idx],
#     all_scores = score_list,
#     best_fit = fit_list[[best_idx]],
#     best_breakpoints = tau_list[[best_idx]],
#     comparison = comparison
#   ))
# }
#
# # Plotting function for the binary segmentation results
# plot_binary_segment_comparison <- function(binary_result) {
#   require(ggplot2)
#   require(dplyr)
#
#   plot_df <- binary_result$comparison_table
#   comparison <- binary_result$comparison_metric
#
#   if (nrow(plot_df) == 0) {
#     stop("No evaluation data to plot")
#   }
#
#   # Determine metric name and direction
#   metric_name <- toupper(comparison)
#   is_lower_better <- comparison %in% c("bic", "loo", "waic", "looic")
#
#   # Find the optimal configuration (matches what the algorithm selected)
#   optimal_breaks <- binary_result$breakpoints
#   optimal_n_breaks <- length(optimal_breaks)
#
#   # Create the plot
#   p <- ggplot(plot_df, aes(x = num_breakpoints, y = metric_value, color = growth_type)) +
#     geom_point(size = 3, alpha = 0.7) +
#     scale_x_continuous(breaks = 0:max(plot_df$num_breakpoints)) +
#     labs(
#       title = paste("Binary Segmentation Model Comparison"),
#       subtitle = paste("Evaluated configurations during recursive search (", metric_name, ")"),
#       x = "Number of Breakpoints",
#       y = metric_name,
#       color = "Growth Type"
#     ) +
#     theme_bw() +
#     theme(
#       legend.position = "bottom",
#       panel.grid.minor = element_blank()
#     )
#
#   # Highlight the selected optimal configuration
#   optimal_row <- plot_df %>%
#     filter(num_breakpoints == optimal_n_breaks) %>%
#     slice_min(metric_value, n = 1)
#
#   if (nrow(optimal_row) > 0) {
#     p <- p +
#       geom_point(data = optimal_row,
#                  aes(x = num_breakpoints, y = metric_value),
#                  color = "black", size = 5, shape = 1, stroke = 2) +
#       annotate("text",
#                x = optimal_row$num_breakpoints[1],
#                y = optimal_row$metric_value[1],
#                label = paste("Selected:", optimal_row$growth_type[1],
#                              "\n", optimal_row$num_breakpoints[1], "breakpoints"),
#                vjust = -1.5, hjust = 0.5, size = 3, fontface = "bold")
#   }
#
#   return(p)
# }
#
#
#
#
#
# plot_trace_parameters <- function(result, model_name, parameters = c("rho", "t0", "K", "n0")) {
#   fit <- result$fits[[model_name]]
#   draws_df <- as_draws_df(fit$draws(variables = parameters))
#   mcmc_trace(draws_df) + ggtitle(paste("Trace Plots -", model_name))
# }
#
# plot_growth_fit <- function(result, model_name, data, breakpoints = numeric(0),
#                             time_grid = NULL, time_col = "time", count_col = "count",
#                             color = "dodgerblue", alpha = 0.3, CI = 0.9) {
#
#   stopifnot(model_name %in% names(result$fits))
#
#   fit <- result$fits[[model_name]]
#   draws <- posterior::as_draws_matrix(fit$draws())
#
#   time_obs <- data[[time_col]]
#   count_obs <- data[[count_col]]
#
#   if (is.null(time_grid)) {
#     time_grid <- seq(min(time_obs), max(time_obs), length.out = 200)
#   }
#
#   G <- length(breakpoints) + 1
#   segment_edges <- c(-Inf, breakpoints, Inf)
#
#   rho_idx <- grep("^rho\\[", colnames(draws))
#   rho_mat <- draws[, rho_idx, drop = FALSE]
#
#   has_t0 <- "t0" %in% colnames(draws)
#   has_n0 <- "n0" %in% colnames(draws)
#
#   if (has_t0) {
#     t0_vec <- draws[, "t0"]
#     n0_vec <- if (has_n0) draws[, "n0"] else rep(1, nrow(draws))
#   } else {
#     t0_vec <- rep(min(time_obs), nrow(draws))
#     n0_vec <- if (has_n0) draws[, "n0"] else rep(1, nrow(draws))
#   }
#
#   K_vec <- if ("K" %in% colnames(draws)) draws[, "K"] else rep(NA, nrow(draws))
#
#   segment_idx <- findInterval(time_grid, vec = segment_edges)
#
#   integrate_r_matrix <- function(t_grid, t0_vec, rho_mat) {
#     n_draws <- nrow(rho_mat)
#     n_times <- length(t_grid)
#     r_mat <- matrix(0, nrow = n_draws, ncol = n_times)
#
#     for (i in seq_along(t_grid)) {
#       t_current <- t_grid[i]
#       for (g in seq_len(G)) {
#         seg_start <- if (g == 1) segment_edges[1] else segment_edges[g]
#         seg_end <- segment_edges[g + 1]
#         if (seg_start >= t_current) next
#
#         if (g == 1) {
#           int_start <- pmax(t0_vec, seg_start)
#         } else {
#           int_start <- rep(seg_start, n_draws)
#         }
#         int_end <- rep(pmin(t_current, seg_end), n_draws)
#
#         valid_idx <- int_start < int_end
#         if (!any(valid_idx)) next
#
#         rho_vec <- rho_mat[, g]
#         if (model_name == "powerlaw") {
#           ratio <- int_end / pmax(int_start, 1e-8)
#           contrib <- rho_vec * log(pmax(ratio, 1e-8))
#         } else {
#           contrib <- rho_vec * (int_end - int_start)
#         }
#
#         contrib[!valid_idx] <- 0
#         r_mat[, i] <- r_mat[, i] + contrib
#       }
#     }
#     return(r_mat)
#   }
#
#   rint_mat <- integrate_r_matrix(time_grid, t0_vec, rho_mat)
#
#   mean_mat <- switch(model_name,
#                      exponential = {
#                        n0_mat <- matrix(n0_vec, nrow = length(n0_vec), ncol = ncol(rint_mat))
#                        n0_mat * exp(rint_mat)
#                      },
#                      powerlaw = {
#                        n0_mat <- matrix(n0_vec, nrow = length(n0_vec), ncol = ncol(rint_mat))
#                        n0_mat * exp(rint_mat)
#                      },
#                      logistic = {
#                        K_mat <- matrix(K_vec, nrow = length(K_vec), ncol = ncol(rint_mat))
#                        n0_mat <- matrix(n0_vec, nrow = length(n0_vec), ncol = ncol(rint_mat))
#
#                        # Ensure K > n0 and no division by ~0
#                        n0_mat <- pmax(n0_mat, 1e-6)
#                        K_mat <- pmax(K_mat, n0_mat + 1e-3)  # enforce K > n0 + ε
#
#                        ratio <- (K_mat - n0_mat) / n0_mat
#                        denom <- 1 + ratio * exp(-rint_mat)
#
#                        denom <- pmax(denom, 1e-6)  # protect against 0
#                        pred <- K_mat / denom
#
#                        # Final cleanup of nonsensical values
#                        pred[!is.finite(pred)] <- NA
#                        pred
#                      },
#                      gompertz = {
#                        K_mat <- matrix(K_vec, nrow = length(K_vec), ncol = ncol(rint_mat))
#                        n0_mat <- matrix(n0_vec, nrow = length(n0_vec), ncol = ncol(rint_mat))
#
#                        loglog_n0_mat <- log(-log(n0_mat / pmax(K_mat, 1e-6)))
#                        exponent_mat <- loglog_n0_mat - rint_mat
#                        K_mat * exp(-exp(exponent_mat))
#                      },
#                      stop("Unknown model: ", model_name)
#   )
#
#   ribbon_df <- apply(mean_mat, 2, function(x) {
#     c(mean = mean(x, na.rm = TRUE),
#       lower = quantile(x, 0.5 - CI / 2, na.rm = TRUE),
#       upper = quantile(x, 0.5 + CI / 2, na.rm = TRUE))
#   }) %>%
#     t() %>%
#     as.data.frame()
#
#   colnames(ribbon_df) <- c("mean", "lower", "upper")
#   ribbon_df$time <- time_grid
#
#   ggplot(ribbon_df, aes(x = time)) +
#     geom_ribbon(aes(ymin = lower, ymax = upper), fill = color, alpha = alpha) +
#     geom_line(aes(y = mean), color = color) +
#     geom_point(data = data, aes_string(x = time_col, y = count_col),
#                color = "black", size = 2) +
#     geom_vline(xintercept = breakpoints, linetype = "dashed", color = "gray") +
#     labs(title = paste("Growth fit:", model_name), x = "Time", y = "Count") +
#     theme_minimal()
# }
#
#
# plot_posterior_growth_rates <- function(result, model_name) {
#   fit <- result$fits[[model_name]]
#   draws <- as_draws_df(fit$draws(variables = "rho"))
#   draws_long <- draws %>%
#     pivot_longer(cols = starts_with("rho["), names_to = "segment", values_to = "value")
#
#   ggplot(draws_long, aes(x = value, fill = segment)) +
#     geom_density(alpha = 0.6) +
#     theme_minimal() +
#     labs(title = paste("Posterior of Growth Rates -", model_name), x = "rho[g]", y = "Density")
# }
#
#
# plot_ppc_ribbons <- function(result, data, time_col = "time", count_col = "count", n_ribbons = 100) {
#   T_vals <- data[[time_col]]
#   N_vals <- data[[count_col]]
#   S <- length(T_vals)
#
#   model_names <- names(result$fits)
#   all_ppc <- lapply(model_names, function(model_name) {
#     fit <- result$fits[[model_name]]
#     yrep <- fit$draws("yrep", format = "draws_matrix")
#
#     summary_df <- apply(yrep, 2, function(x) {
#       c(mean = mean(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975))
#     })
#     rownames(summary_df) = c("mean", "lower", "upper")
#
#     df <- data.frame(t(summary_df))
#     df$time <- T_vals
#     df$obs <- N_vals
#     df$model <- model_name
#     df
#   })
#
#   ppc_df <- do.call(rbind, all_ppc)
#
#   ggplot(ppc_df, aes(x = time)) +
#     geom_ribbon(aes(ymin = lower, ymax = upper), fill = "skyblue", alpha = 0.2) +
#     geom_line(aes(y = mean), color = "blue") +
#     geom_point(aes(y = obs), color = "black", size = 1) +
#     facet_wrap(~ model, scales = "free_y") +
#     labs(title = "Posterior Predictive Checks for All Models",
#          y = "Count", x = "Time") +
#     theme_bw()
# }
#
# plot_posterior_t0 <- function(result, model_name, overlay_fit = FALSE,
#                               data = NULL, breakpoints = numeric(0)) {
#   fit <- result$fits[[model_name]]
#   if (!"t0" %in% colnames(as_draws_matrix(fit$draws()))) {
#     warning("This model does not estimate t0.")
#     return(NULL)
#   }
#   t0_vals <- posterior::as_draws_df(fit$draws("t0"))$t0
#
#   p <- ggplot(data.frame(t0 = t0_vals), aes(x = t0)) +
#     geom_density(fill = "salmon", alpha = 0.6) +
#     labs(title = paste("Posterior of Initiation Time (t0) -", model_name),
#          x = "t0", y = "Density") +
#     theme_minimal()
#
#   if (overlay_fit && !is.null(data)) {
#     fit_plot <- plot_growth_fit(result, model_name, data, breakpoints = breakpoints)
#     p <- cowplot::plot_grid(fit_plot, p, ncol = 1, rel_heights = c(2, 1))
#   }
#   return(p)
# }
#
# plot_best_model_ribbon <- function(result, data, time_col = "time", count_col = "count") {
#   if (!("model_table" %in% names(result) && "fits" %in% names(result))) {
#     stop("Result must include 'model_table' and 'fits'.")
#   }
#
#   metric_col <- colnames(result$model_table)[1]
#   metric <- tolower(metric_col)
#   is_logml <- grepl("log", metric)
#
#   # Pick best model: min for BIC/LOO/WAIC, max for log marginal likelihood
#   best_model <- rownames(result$model_table)[which(
#     if (is_logml) result$model_table[[1]] == max(result$model_table[[1]])
#     else result$model_table[[1]] == min(result$model_table[[1]])
#   )][1]
#
#   fit <- result$fits[[best_model]]
#
#   T_vals <- data[[time_col]]
#   N_vals <- data[[count_col]]
#   yrep <- fit$draws("yrep", format = "draws_matrix")
#
#   summary_df <- apply(yrep, 2, function(x) {
#     c(mean = mean(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975))
#   })
#   rownames(summary_df) = c("mean", "lower", "upper")
#   df <- as.data.frame(t(summary_df))
#   df$time <- T_vals
#   df$obs <- N_vals
#
#   ggplot(df, aes(x = time)) +
#     geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightgreen", alpha = 0.2) +
#     geom_line(aes(y = mean), color = "forestgreen") +
#     geom_point(aes(y = obs)) +
#     labs(
#       title = paste("Best Model:", best_model, "(", toupper(metric), ")"),
#       y = "Count", x = "Time"
#     ) +
#     theme_minimal()
# }
#
# plot_model_comparison <- function(result) {
#   df <- result$model_table
#   metric_col <- colnames(df)[1]
#   metric_name <- metric_col
#   is_logml <- grepl("log", metric_name, ignore.case = TRUE)
#
#   df$model <- rownames(df)
#
#   # Sort models: lower is better for most (BIC, LOO, WAIC), higher is better for log marginal likelihood
#   df <- df[order(if (is_logml) -df[[metric_col]] else df[[metric_col]]), ]
#   df$model <- factor(df$model, levels = df$model)  # preserve sort order for plotting
#
#   ggplot(df, aes(x = model, y = .data[[metric_col]])) +
#     geom_point(size = 3) +
#     geom_text(aes(label = round(.data[[metric_col]], 1)), vjust = -0.5, size = 3.5) +
#     labs(title = paste("Model Comparison (", toupper(metric_name), ")", sep = ""),
#          x = "Model", y = metric_name) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
# }
#
# fit_best_recovery_model <- function(data,
#                                     stan_dir,
#                                     chains = 4,
#                                     iter = 4000,
#                                     seed = 123,
#                                     cores = 4,
#                                     comparison = c("loo", "bic")) {
#   comparison <- match.arg(comparison)
#
#   stopifnot(all(c("time", "count") %in% colnames(data)))
#   data <- data[order(data$time), ]
#
#   stan_data <- list(
#     S = nrow(data),
#     N = data$count,
#     T = data$time
#   )
#
#   # Model files
#   model_files <- c(
#     "two_pop_denovo.stan",
#     "two_pop_preexisting.stan",
#     "two_pop_single.stan",
#     "two_pop_both.stan"
#   )
#
#   fits <- list()
#   ic_values <- numeric(length(model_files))
#
#   for (i in seq_along(model_files)) {
#     mod <- cmdstan_model(file.path(stan_dir, model_files[i]))
#     fit <- mod$sample(
#       data = stan_data,
#       chains = chains,
#       iter_warmup = iter,
#       iter_sampling = iter,
#       seed = seed,
#       parallel_chains = cores
#     )
#     fits[[i]] <- fit
#
#     # Compute IC
#     log_lik <- fit$draws("log_lik") |>
#       posterior::as_draws_matrix()
#     if (comparison == "loo") {
#       ic_values[i] <- loo(log_lik)$estimates["elpd_loo", "Estimate"]
#     } else {
#       # Extract log posterior likelihood
#       log_lik_sum <- sum(apply(log_lik, 1, mean))
#       k <- length(fit$metadata()$parameters) # approx param count
#       n <- nrow(data)
#       ic_values[i] <- -2 * log_lik_sum + k * log(n)  # BIC
#     }
#   }
#
#   # Select best
#   if (comparison == "loo") {
#     best_idx <- which.max(ic_values) # higher elpd is better
#   } else {
#     best_idx <- which.min(ic_values) # lower BIC is better
#   }
#
#   list(
#     best_model = model_files[best_idx],
#     best_fit = fits[[best_idx]],
#     all_fits = fits,
#     scores = ic_values
#   )
# }
#
# plot_recovery_model_fit <- function(fit_result, data, model_label = NULL,
#                                     time_col = "time", count_col = "count",
#                                     CI = 0.9, color = "dodgerblue", alpha = 0.3) {
#   stopifnot(all(c(time_col, count_col) %in% colnames(data)))
#
#   fit <- fit_result$best_fit
#   if (is.null(model_label)) model_label <- fit_result$best_model
#
#   # Extract draws for yrep (posterior predictive means)
#   draws <- as_draws_matrix(fit$draws("yrep"))
#   time_obs <- data[[time_col]]
#   count_obs <- data[[count_col]]
#
#   # Compute credible intervals
#   pred_summary <- apply(draws, 2, quantile, probs = c((1 - CI) / 2, 0.5, 1 - (1 - CI) / 2))
#   pred_df <- data.frame(
#     time = time_obs,
#     lower = pred_summary[1, ],
#     median = pred_summary[2, ],
#     upper = pred_summary[3, ]
#   )
#
#   # Plot
#   ggplot() +
#     geom_point(aes(x = time_obs, y = count_obs), color = "black", size = 2) +
#     geom_ribbon(aes(x = pred_df$time, ymin = pred_df$lower, ymax = pred_df$upper),
#                 fill = color, alpha = alpha) +
#     geom_line(aes(x = pred_df$time, y = pred_df$median), color = color, size = 1.2) +
#     labs(title = paste("Best Fitted Model:", model_label),
#          x = "Time",
#          y = "Observed & Predicted Counts") +
#     theme_minimal(base_size = 14)
# }
#
#
#
# detect_breakpoints_rss <- function(data, K, min_segment_size = NULL) {
#   df <- data.frame(time = data$time, y = log1p(data$count))
#   N <- nrow(df)
#
#   if (is.null(min_segment_size)) {
#     min_segment_size <- max(5, floor(N / (K * 3)))  # Ensure reasonable segment sizes
#   }
#
#   # Function to calculate RSS for a given set of breakpoints
#   calculate_rss <- function(breakpoints) {
#     bp_full <- c(1, breakpoints, N + 1)
#     total_rss <- 0
#
#     for (i in 1:(length(bp_full) - 1)) {
#       start_idx <- bp_full[i]
#       end_idx <- bp_full[i + 1] - 1
#
#       if (end_idx - start_idx + 1 >= 3) {  # Need at least 3 points for regression
#         segment_data <- df[start_idx:end_idx, ]
#         fit <- lm(y ~ time, data = segment_data)
#         total_rss <- total_rss + sum(residuals(fit)^2)
#       } else {
#         total_rss <- total_rss + 1e6  # Penalty for too-small segments
#       }
#     }
#     return(total_rss)
#   }
#
#   # Dynamic programming approach for optimal breakpoints
#   if (K == 1) {
#     return(numeric(0))
#   }
#
#   # Generate candidate breakpoints (exclude too-close-to-edges points)
#   candidates <- seq(min_segment_size + 1, N - min_segment_size, by = 1)
#
#   # Use dynamic programming or greedy search for optimal breakpoints
#   best_breakpoints <- numeric(0)
#   remaining_candidates <- candidates
#
#   for (k in 1:(K - 1)) {
#     best_rss <- Inf
#     best_bp <- NULL
#
#     for (candidate in remaining_candidates) {
#       test_bp <- sort(c(best_breakpoints, candidate))
#
#       # Check minimum separation constraint
#       if (length(test_bp) == 1 || all(diff(c(1, test_bp, N + 1)) >= min_segment_size)) {
#         rss <- calculate_rss(test_bp)
#         if (rss < best_rss) {
#           best_rss <- rss
#           best_bp <- candidate
#         }
#       }
#     }
#
#     if (!is.null(best_bp)) {
#       best_breakpoints <- c(best_breakpoints, best_bp)
#       # Remove candidates too close to selected breakpoint
#       remaining_candidates <- remaining_candidates[abs(remaining_candidates - best_bp) >= min_segment_size]
#     }
#   }
#
#   return(sort(best_breakpoints))
# }
#
# fit_lr_hmm_improved <- function(data, K = 3, max_iter = 100, tol = 1e-10,
#                                 verbose = F, init_method = "rss", min_segment_size = 5,
#                                 blend_width = 0.15, temporal_smooth = 0.05) {
#
#   logSumExp <- function(x) {
#     m <- max(x); if (!is.finite(m)) return(-Inf); m + log(sum(exp(x - m)))
#   }
#
#   stopifnot(all(c("time", "count") %in% colnames(data)))
#   df <- data.frame(time = data$time, y = log1p(data$count))
#   N <- nrow(df)
#
#   ## --- STEP 1: Improved breakpoint initialization ---
#   bp_idx = detect_breakpoints_rss(data, K, min_segment_size = min_segment_size)
#
#   # Ensure we have the right number of breakpoints
#   if (length(bp_idx) < K - 1) {
#     # Fallback: supplement with quantile-based breakpoints
#     missing <- (K - 1) - length(bp_idx)
#     quantile_bp <- floor(quantile(1:N, probs = seq(0, 1, length.out = K + missing + 1)))[2:(K + missing)]
#     # Remove any that are too close to existing breakpoints
#     min_sep <- floor(N / (K * 2))
#     for (qbp in quantile_bp) {
#       if (length(bp_idx) == 0 || all(abs(bp_idx - qbp) >= min_sep)) {
#         bp_idx <- c(bp_idx, qbp)
#         if (length(bp_idx) >= K - 1) break
#       }
#     }
#   } else if (length(bp_idx) > K - 1) {
#     # Too many breakpoints, keep the best K-1
#     bp_idx <- bp_idx[1:(K-1)]
#   }
#
#   bp_idx <- sort(unique(pmax(2, pmin(N - 1, bp_idx))))
#   bp_times <- df$time[bp_idx]
#
#   if (verbose) {
#     cat("Initial breakpoints at indices:", paste(bp_idx, collapse = ", "), "\n")
#     cat("Initial breakpoints at times:", paste(round(bp_times, 3), collapse = ", "), "\n")
#   }
#
#   ## --- STEP 2: HMM transition matrix (adaptive to breakpoint quality) ---
#   if (K == 1) {
#     A <- matrix(1, 1, 1); pi <- 1
#   } else {
#     seg_lengths <- diff(c(1, bp_idx, N))
#     mean_seg_length <- mean(seg_lengths)
#
#     # Adaptive stay probabilities based on segment characteristics
#     stay_prob <- pmin(0.95, pmax(0.6, 1 - 2 / (1 + seg_lengths / mean_seg_length)))
#
#     A <- matrix(0, K, K)
#     diag(A) <- stay_prob
#
#     # Allow transitions to adjacent states
#     for (k in 1:(K-1)) {
#       A[k, k+1] <- 1 - stay_prob[k]
#     }
#     A[K, K] <- 1  # Last state is absorbing
#
#     # Initial state probabilities favor starting in first state
#     pi <- numeric(K)
#     pi[1] <- 0.8
#     if (K > 1) pi[2:K] <- 0.2 / (K - 1)
#     pi <- pi / sum(pi)
#   }
#
#   log_pi <- log(pi)
#   logA <- matrix(-Inf, K, K)
#   logA[A > 0] <- log(A[A > 0])
#
#   ## --- STEP 3: Smart responsibility initialization ---
#   responsibilities <- matrix(0, N, K)
#
#   if (K == 1) {
#     responsibilities[, 1] <- 1
#   } else {
#     # Create hard assignments based on breakpoints
#     seg_bounds <- c(1, bp_idx, N)
#
#     for (k in 1:K) {
#       start_idx <- seg_bounds[k]
#       end_idx <- seg_bounds[k + 1]
#       responsibilities[start_idx:end_idx, k] <- 1
#     }
#
#     # Apply soft blending near breakpoints (adaptive width)
#     blend_points <- floor(blend_width * N / K)  # Blend width as fraction of typical segment
#
#     for (i in 1:length(bp_idx)) {
#       bp <- bp_idx[i]
#       width <- min(blend_points, floor(min(seg_lengths) / 3))  # Don't exceed 1/3 of smallest segment
#
#       if (width > 0) {
#         blend_range <- max(1, bp - width):min(N, bp + width)
#
#         for (t in blend_range) {
#           if (t <= N) {
#             # Distance-based blending
#             dist_factor <- 1 - abs(t - bp) / width
#             blend_strength <- plogis(3 * (dist_factor - 0.5))  # Sigmoid blending
#
#             # Redistribute probability between adjacent states
#             if (i < K) {  # Not the last breakpoint
#               orig_left <- responsibilities[t, i]
#               orig_right <- responsibilities[t, i + 1]
#
#               responsibilities[t, i] <- orig_left * (1 - 0.4 * blend_strength) +
#                 orig_right * (0.4 * blend_strength)
#               responsibilities[t, i + 1] <- orig_right * (1 - 0.4 * blend_strength) +
#                 orig_left * (0.4 * blend_strength)
#             }
#           }
#         }
#       }
#     }
#   }
#
#   # Normalize responsibilities
#   responsibilities <- responsibilities / rowSums(responsibilities)
#
#   ## --- STEP 4: Temporal smoothing (adaptive) ---
#   if (temporal_smooth > 0) {
#     smooth_w <- max(2, floor(temporal_smooth * N))
#     temp_resp <- responsibilities
#
#     for (t in (smooth_w + 1):(N - smooth_w)) {
#       window_resp <- responsibilities[(t - smooth_w):(t + smooth_w), , drop = FALSE]
#       # Weighted average with center emphasis
#       weights <- dnorm(-smooth_w:smooth_w, 0, smooth_w/2)
#       weights <- weights / sum(weights)
#
#       temp_resp[t, ] <- apply(window_resp, 2, function(col) sum(col * weights))
#     }
#
#     # Blend smoothed with original (preserve breakpoint structure)
#     alpha <- 0.3  # Smoothing strength
#     responsibilities <- (1 - alpha) * responsibilities + alpha * temp_resp
#     responsibilities <- responsibilities / rowSums(responsibilities)
#   }
#
#   ## --- STEP 5: EM Algorithm with monitoring ---
#   prev_loglik <- -Inf
#   loglik_history <- numeric(max_iter)
#
#   for (iter in 1:max_iter) {
#     if (verbose && iter %% 10 == 1) cat("Iteration", iter, "\n")
#
#     ## M-step: Fit linear models for each state
#     models <- lapply(1:K, function(k) {
#       w <- responsibilities[, k]
#       w_sum <- sum(w)
#
#       if (w_sum < 1e-3) {
#         # Degenerate state - use global fit as fallback
#         global_fit <- lm(y ~ time, df)
#         return(list(
#           alpha = coef(global_fit)[1],
#           beta = coef(global_fit)[2],
#           sigma = summary(global_fit)$sigma
#         ))
#       }
#
#       # Weighted regression
#       fit <- tryCatch({
#         lm(y ~ time, df, weights = w)
#       }, error = function(e) {
#         lm(y ~ time, df)  # Fallback to unweighted
#       })
#
#       if (is.null(fit) || any(is.na(coef(fit)))) {
#         return(list(alpha = 0, beta = 0, sigma = 1))
#       }
#
#       # Robust sigma estimation
#       residuals_weighted <- residuals(fit) * sqrt(w)
#       sigma <- sqrt(max(sum(residuals_weighted^2) / max(w_sum - 2, 1), 1e-4))
#
#       list(
#         alpha = coef(fit)[1],
#         beta = coef(fit)[2],
#         sigma = sigma
#       )
#     })
#
#     ## E-step: Forward-backward algorithm
#     log_em <- matrix(-Inf, N, K)
#     for (k in 1:K) {
#       mu <- models[[k]]$alpha + models[[k]]$beta * df$time
#       s <- models[[k]]$sigma
#       if (s > 1e-6) {
#         log_em[, k] <- dnorm(df$y, mu, s, log = TRUE)
#       }
#     }
#
#     # Forward pass
#     log_alpha <- matrix(-Inf, N, K)
#     log_alpha[1, ] <- log_pi + log_em[1, ]
#
#     for (t in 2:N) {
#       for (k in 1:K) {
#         # Only consider valid previous states
#         valid_prev <- which(logA[, k] > -Inf)
#         if (length(valid_prev) > 0) {
#           log_alpha[t, k] <- log_em[t, k] +
#             logSumExp(log_alpha[t - 1, valid_prev] + logA[valid_prev, k])
#         }
#       }
#     }
#
#     # Backward pass
#     log_beta <- matrix(-Inf, N, K)
#     log_beta[N, ] <- 0
#
#     for (t in (N - 1):1) {
#       for (k in 1:K) {
#         valid_next <- which(logA[k, ] > -Inf)
#         if (length(valid_next) > 0) {
#           log_beta[t, k] <- logSumExp(logA[k, valid_next] +
#                                         log_em[t + 1, valid_next] +
#                                         log_beta[t + 1, valid_next])
#         }
#       }
#     }
#
#     # Update responsibilities
#     log_gamma <- log_alpha + log_beta
#     responsibilities <- exp(sweep(log_gamma, 1, apply(log_gamma, 1, logSumExp), "-"))
#
#     # Handle numerical issues
#     responsibilities[is.nan(responsibilities)] <- 1/K
#     responsibilities <- responsibilities / rowSums(responsibilities)
#
#     # Calculate log-likelihood
#     loglik <- logSumExp(log_alpha[N, ])
#     loglik_history[iter] <- loglik
#
#     if (verbose && iter %% 10 == 1) cat("Log-likelihood:", round(loglik, 4), "\n")
#
#     # Check convergence
#     if (iter > 1 && abs(loglik - prev_loglik) < tol) {
#       if (verbose) cat("Converged at iteration", iter, "\n")
#       break
#     }
#
#     # Check for likelihood decrease (potential numerical issues)
#     if (iter > 5 && loglik < prev_loglik - abs(prev_loglik) * 1e-6) {
#       if (verbose) cat("Warning: Likelihood decreased at iteration", iter, "\n")
#     }
#
#     prev_loglik <- loglik
#   }
#
#   ## --- STEP 6: Final state assignment and breakpoint extraction ---
#   final_states <- apply(responsibilities, 1, which.max)
#
#   # Extract refined breakpoints from state transitions
#   state_changes <- which(diff(final_states) != 0)
#   refined_breakpoints <- df$time[state_changes]
#
#   # Calculate growth rates and other statistics
#   growth_rates <- sapply(models, function(m) m$beta)
#
#   # Model diagnostics
#   state_proportions <- colMeans(responsibilities)
#
#   list(
#     data = df,
#     states = final_states,
#     responsibilities = responsibilities,
#     models = models,
#     loglik = loglik,
#     loglik_history = loglik_history[1:iter],
#     transition_matrix = A,
#     initial_probabilities = pi,
#     initial_breakpoints = bp_times,
#     refined_breakpoints = refined_breakpoints,
#     growth_rates = growth_rates,
#     state_proportions = state_proportions,
#     convergence_iter = iter,
#     initialization_method = init_method
#   )
# }
#
# segment_fit <- function(data,
#                         max_segments = 4,
#                         min_segment_size = 3,
#                         comparison = "bic",
#                         enforce_rho_separation = TRUE,
#                         alpha_rho = 0.05,
#                         models_to_fit = c("exponential"),
#                         ...) {
#   require(dplyr)
#   require(posterior)
#
#   rhos_distinct <- function(draws, alpha = 0.05) {
#     rho_idx <- grep("^rho\\[", colnames(draws))
#     rho_mat <- draws[, rho_idx, drop = FALSE]
#     G <- ncol(rho_mat)
#     if (G <= 1) return(TRUE)
#     ci_bounds <- apply(rho_mat, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2))
#     for (g in seq_len(G - 1)) {
#       if (!(ci_bounds[2, g] < ci_bounds[1, g + 1] || ci_bounds[2, g + 1] < ci_bounds[1, g])) {
#         return(FALSE)
#       }
#     }
#     TRUE
#   }
#
#   # --- wrapper for growth model evaluation ---
#   evaluate_model <- function(breaks) {
#     fit_growth_models(
#       data = data,
#       breakpoints = breaks,
#       comparison = comparison,
#       models_to_fit = models_to_fit,
#       ...
#     )
#   }
#
#   # evaluate_model <- function(breaks) {
#   #   fit_growth_models(
#   #     data = data,
#   #     breakpoints = breaks,
#   #     comparison = comparison,
#   #     models_to_fit = models_to_fit, with_initiation = with_initiation, stan_dir = stan_dir, chains = chains, iter = iter, seed = seed, cores = cores
#   #   )
#   # }
#
#   store_evaluation <- function(breaks, result) {
#     metric_col <- if (comparison == "bic") "BIC" else "looic"
#     for (growth_type in rownames(result$model_table)) {
#       all_evaluations <<- append(all_evaluations, list(data.frame(
#         num_breakpoints = length(breaks),
#         num_segments = length(breaks) + 1,
#         breakpoints = if (length(breaks) == 0) "none" else paste(round(breaks, 3), collapse = ", "),
#         growth_type = growth_type,
#         metric_value = result$model_table[growth_type, metric_col],
#         stringsAsFactors = FALSE
#       )))
#     }
#   }
#
#   all_evaluations <- list()
#   all_checked_candidates <- list()
#
#   results_list <- lapply(1:max_segments, function(k) {
#     hmm_fit = fit_lr_hmm_improved(data = data, K = k, min_segment_size = min_segment_size)
#     proposed_bp = hmm_fit$refined_breakpoints
#     #proposed_bp <- data$time[which(diff(hmm_fit$states) != 0)]
#
#     print(proposed_bp)
#
#     # Check segment size constraint
#     segments <- sort(c(min(data$time), proposed_bp, max(data$time)))
#     segment_lengths <- sapply(seq_along(segments[-1]), function(i) {
#       sum(data$time >= segments[i] & data$time < segments[i + 1])
#     })
#
#     valid <- all(segment_lengths >= min_segment_size)
#     rejected_reason <- if (!valid) "min_segment_size" else NA
#
#     if (valid) {
#       result <- evaluate_model(proposed_bp)
#       store_evaluation(proposed_bp, result)
#
#       if (enforce_rho_separation) {
#         model_draws <- posterior::as_draws_matrix(result$fits[[1]]$draws())
#         if (!rhos_distinct(model_draws, alpha = alpha_rho)) {
#           return(list(breakpoints = proposed_bp, num_segments = k,
#                       valid = FALSE, rejected_reason = "rho_overlap",
#                       result = result, metric = NA))
#         }
#       }
#
#       metric <- if (comparison == "bic") min(result$model_table$BIC) else result$model_table[1, "looic"]
#       return(list(breakpoints = proposed_bp, num_segments = k,
#                   valid = TRUE, rejected_reason = NA,
#                   result = result, metric = metric))
#     } else {
#       return(list(breakpoints = proposed_bp, num_segments = k,
#                   valid = FALSE, rejected_reason = rejected_reason,
#                   result = NULL, metric = NA))
#     }
#   })
#
#   all_checked_candidates <- results_list
#   valid_evals <- Filter(function(x) isTRUE(x$valid) && !is.na(x$metric), results_list)
#
#   if (length(valid_evals) == 0) stop("No valid candidate models were found.")
#
#   best_candidate <- valid_evals[[which.min(sapply(valid_evals, function(x) x$metric))]]
#   best_candidate$comparison_table <- do.call(rbind, all_evaluations)
#   best_candidate$comparison_metric <- comparison
#   best_candidate$all_checked_candidates <- all_checked_candidates
#
#   return(best_candidate)
# }
#
