#' #' Fit growth model to bipod object
#' #'
#' #' @param x a bipod object
#' #' @param factor_size numeric factor by which to divide counts in the bipod object
#' #' @param available_changepoints .
#' #' @param max_iter .
#' #' @param min_support_points .
#' #'
#' #' @return the input bipod object with an added 'breakpoints_fit' slot containing the fitted model for the breakpoints
#' #' @export
#' breakpoints_inference <- function(
#'     x,
#'     factor_size = 1,
#'     available_changepoints = c(0:2),
#'     min_support_points = 2,
#'     max_iter = 20) {
#'   # Check input
#'   if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
#'   if (!(factor_size > 0)) stop("factor_size must be positive")
#'
#'   # Get the model
#'   model_name <- "piecewise_changepoints"
#'   m <- biPOD:::get_model(model_name = model_name)
#'
#'   # Loop
#'   res <- dplyr::tibble()
#'
#'   fits <- list()
#'   cli::cli_progress_bar("Inferring changepoints inference", total = max_iter * length(available_changepoints), clear = F)
#'   for (J in available_changepoints) {
#'
#'     if (J >= 1) {
#'       input_data <- list(
#'         S = nrow(x$counts),
#'         G = J,
#'         N = x$counts$count / factor_size,
#'         T = x$counts$time - min(x$counts$time),
#'         b_prior = array(biPOD:::find_equispaced_points(min(x$counts$time - min(x$counts$time)), max(x$counts$time - min(x$counts$time)), N = J), dim = c(J)),
#'         sigma_changepoints = max(x$counts$time) - min(x$counts$time)
#'       )
#'     } else {
#'       input_data <- list(
#'         S = nrow(x$counts),
#'         G = J,
#'         N = x$counts$count / factor_size,
#'         T = x$counts$time - min(x$counts$time),
#'         b_prior = array(0, dim = c(0)),
#'         sigma_changepoints = max(x$counts$time) - min(x$counts$time)
#'       )
#'     }
#'
#'     for (i in 1:max_iter) {
#'       lp <- -Inf
#'
#'       #out <- utils::capture.output(suppressMessages(f_pf <- m$pathfinder(input_data)))
#'       #f_pf$output()
#'
#'       out <- tryCatch({
#'         out <- utils::capture.output(suppressMessages(f_pf <- m$pathfinder(input_data)))
#'         k <- ncol(f_pf$draws()) - 3 # number of parameters
#'         n <- nrow(x$counts)         # number of obs
#'         lp <- f_pf$draws("lp__") %>% stats::median()
#'         bic <- k * log(n) - 2 * lp
#'         f_pf
#'       }, error = function(cond) {
#'         lp <- Inf
#'         bic <- Inf
#'         return(NA)
#'       }, warning = function(warn) {
#'         lp <- Inf
#'         bic <- Inf
#'         return(NA)
#'       }
#'       )
#'
#'       res <- dplyr::bind_rows(res, dplyr::tibble(J = J, iter = i, lp = lp, bic = bic))
#'       fits[[paste0(J, " _ ", i)]] <- out
#'
#'       cli::cli_progress_update()
#'     }
#'   }
#'
#'   res <- res %>%
#'     dplyr::filter(bic != Inf) %>%
#'     stats::na.omit()
#'
#'   status <- 'FAIL'
#'   while ((status == 'FAIL') & (nrow(res) > 0)) {
#'     best_res <- res %>% dplyr::filter(bic == min(bic))
#'
#'     best_fit <- fits[[paste0(best_res$J, " _ ", best_res$iter)]]
#'
#'     if (best_res$J == 0) {
#'       status <- 'PASS'
#'     } else {
#'       draws <- best_fit$draws(format = "matrix")
#'       b_draws <- draws[,grepl("b", colnames(draws))] %>% dplyr::as_tibble()
#'
#'       median_breakpoints <- b_draws %>%
#'         dplyr::summarise_all(stats::median) %>%
#'         as.numeric()
#'
#'       low_bp <- b_draws %>%
#'         dplyr::summarise_all(low_quant) %>%
#'         as.numeric()
#'
#'       high_bp <- b_draws %>%
#'         dplyr::summarise_all(high_quant) %>%
#'         as.numeric()
#'
#'       median_breakpoints <- median_breakpoints + min(x$counts$time)
#'       low_bp <- low_bp + min(x$counts$time)
#'       high_bp <- high_bp + min(x$counts$time)
#'
#'       groups <- bp_to_groups(x$counts, median_breakpoints)
#'       min_group_numerosity <- dplyr::tibble(groups = groups) %>%
#'         dplyr::group_by(groups) %>%
#'         dplyr::summarise(n = dplyr::n()) %>%
#'         dplyr::pull(.data$n) %>%
#'         min()
#'
#'       if (all(high_bp <= max(x$counts$time)) & all(low_bp >= min(x$counts$time)) & (min_group_numerosity > min_support_points)) {
#'         status <- 'PASS'
#'       } else {
#'         res <- res %>% dplyr::filter(lp != max(lp))
#'       }
#'     }
#'   }
#'
#'   # Prep final fit
#'   if (best_res$J >= 1) {
#'     draws <- best_fit$draws(format = "matrix")
#'     b_draws <- draws[,grepl("b", colnames(draws))] %>% dplyr::as_tibble()
#'
#'     median_breakpoints <- b_draws %>%
#'       dplyr::summarise_all(stats::median) %>%
#'       as.numeric()
#'
#'     sd_breakpoints <- b_draws %>%
#'       dplyr::summarise_all(stats::sd) %>%
#'       as.numeric()
#'
#'     input_data <- list(
#'       S = nrow(x$counts),
#'       G = best_res$J,
#'       N = x$counts$count / factor_size,
#'       T = x$counts$time - min(x$counts$time),
#'       b_prior = median_breakpoints,
#'       sigma_changepoints = max(sd_breakpoints)
#'     )
#'   } else {
#'     input_data <- list(
#'       S = nrow(x$counts),
#'       G = best_res$J,
#'       N = x$counts$count / factor_size,
#'       T = x$counts$time - min(x$counts$time),
#'       b_prior = array(0, dim = c(0)),
#'       sigma_changepoints = .01
#'     )
#'   }
#'   tmp <- utils::capture.output(suppressMessages(best_fit <- final_fit <- m$sample(data=input_data, chains=4, parallel_chains=4, iter_warmup = 4000, iter_sampling = 4000)))
#'
#'   # Store results
#'   elbo_data <- c()
#'   # if (variational) elbo_data <- elbo_d %>% stats::na.omit()
#'   # fit <- fit_model
#'
#'   # Add results to bipod object
#'   x$breakpoints_elbo <- elbo_data
#'   x$breakpoints_fit <- convert_mcmc_fit_to_biPOD(best_fit)
#'
#'   # Write fit info
#'   # x$metadata$sampling <- sampling
#'   x$metadata$factor_size <- factor_size
#'   # x$metadata$prior_K <- input_data$prior_K
#'
#'   # Add median of breakpoints
#'   # n_changepoints <- length(input_data$changing_times_prior)
#'   # breakpoints_names <- lapply(1:n_changepoints, function(i) {
#'   #   paste0("changing_times[", i, "]")
#'   # }) %>% unlist()
#'
#'   if (best_res$J == 0) {
#'     median_breakpoints = NULL
#'   } else {
#'     median_breakpoints <- best_fit$draws(variables = 'b', format = 'matrix') %>%
#'       dplyr::as_tibble() %>%
#'       dplyr::summarise_all(stats::median) %>%
#'       as.numeric()
#'
#'     median_breakpoints <- median_breakpoints + min(x$counts$time)
#'   }
#'
#'   x$metadata$breakpoints <- median_breakpoints
#'
#'   if (!(is.null(median_breakpoints))) {
#'     x$counts$group <- bp_to_groups(x$counts, x$metadata$breakpoints)
#'   }
#'
#'   cli::cli_alert_success("Breakpoints have been inferred. Inspect the results using the {.field plot_breakpoints_posterior} function.")
#'   cli::cli_alert_info("Median of the inferred breakpoints have been succesfully stored.")
#'
#'   x
#' }
#'
#'
#' ## Utils for breakpoints inference
#'
#' # prep_data_bp_inference <- function(x, factor_size, G) {
#' #
#' #   if (G >= 1) {
#' #     input_data <- list(
#' #       S = nrow(x$counts),
#' #       G = G,
#' #       N = as.array(x$counts$count / factor_size),
#' #       T = as.array(x$counts$time - min(x$counts$time)),
#' #       b_prior = find_equispaced_points(min(d$time - min(d$time)), max(d$time - min(d$time)), N = G),
#' #       sigma_changepoints = 10
#' #     )
#' #   } else {
#' #     input_data <- list(
#' #       S = nrow(x$counts),
#' #       G = G,
#' #       N = as.array(x$counts$count / factor_size),
#' #       T = as.array(x$counts$time - min(x$counts$time)),
#' #       b_prior = array(0, dim = c(0)),
#' #       sigma_changepoints = 10
#' #     )
#' #   }
#' #   return(input_data)
#' # }
#'
#'
#' find_equispaced_points <- function(start, end, N) {
#'   if (N < 1) {
#'     stop("N must be greater than 1")
#'   }
#'
#'   equispaced_points <- seq(start, end, length.out = N + 2)[-c(1, N + 2)]
#'   return(equispaced_points)
#' }
#'
#' low_quant = function(x) { stats::quantile(x, 0.05) }
#' high_quant = function(x) { stats::quantile(x, 0.95) }
#'
#'
#' run_pathfinder = function(data,mod){
#'   data_file = tempfile(fileext='.json')
#'   output_file = tempfile(fileext='.csv')
#'   cmdstanr::write_stan_json(data,file=data_file)
#'   processx::run(
#'     command = mod$exe_file()
#'     , args = c(
#'       'pathfinder'
#'       , 'data'
#'       , paste0('file=',data_file)
#'       , 'output'
#'       , paste0('file=',output_file)
#'     )
#'     , echo_cmd = T
#'     , stdout = ""
#'     , stderr = "2>&1"
#'   )
#'
#'   (
#'     paste0("grep '^[#l]' '",output_file,"'")
#'     %>% system(intern=T)
#'     %>% strsplit('\n')
#'     %>% unlist()
#'   ) -> header
#'   header_nlines = which(stringr::str_starts(header,'lp'))
#'   found_samples_col_names = unlist(strsplit(header[header_nlines],','))
#'   (
#'     data.table::fread(
#'       cmd = paste0(
#'         "tail -n+"
#'         , header_nlines + 1
#'         , " '"
#'         , output_file
#'         , "' | grep -v '^[#]' --color=never"
#'       )
#'       , data.table = FALSE
#'       , sep = ','
#'       , header = F
#'       , col.names = found_samples_col_names
#'       , colClasses = list(numeric=1:length(found_samples_col_names))
#'     )
#'     %>% dplyr::as_tibble()
#'   ) ->
#'     out
#'   return(out)
#' }
