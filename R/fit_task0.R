#' Fit growth model to bipod object
#'
#' @param x a bipod object
#' @param norm .
#' @param n_trials .
#' @param avg_points_per_window .
#' @param max_breakpoints .
#'
#' @return the input bipod object with an added 'breakpoints_fit' slot containing the fitted model for the breakpoints
#' @export
fit_breakpoints = function(
  x,
  norm=TRUE,
  n_trials=500,
  avg_points_per_window = 3,
  max_breakpoints = 10
  ) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  #if (!(factor_size > 0)) stop("factor_size must be positive")
  d <- x$counts

  res <- find_breakpoints(
    d,
    avg_points_per_window=avg_points_per_window,
    max_breakpoints=max_breakpoints,
    norm=norm,
    n_trials=n_trials
  )

  best_bp <- res$best_bp
  best_fit <- res$best_fit

  # Store results
  elbo_data <- c()
  # if (variational) elbo_data <- elbo_d %>% stats::na.omit()
  # fit <- fit_model

  # Add results to bipod object
  #x$breakpoints_elbo <- elbo_data
  x$breakpoints_fit <- best_fit

  x$metadata$breakpoints <- best_bp

  if (!(is.null(best_bp))) {
    x$counts$group <- bp_to_groups(x$counts, x$metadata$breakpoints)
  }
  cli::cli_alert_info("Median of the inferred breakpoints have been succesfully stored.")

  x

}


find_breakpoints = function(d, avg_points_per_window, max_breakpoints, norm, n_trials) {
  x <- d$time
  y <- log(d$count)

  if (norm) {
    x <- (x - mean(x)) / stats::sd(x)
    y <- (y - mean(y)) / stats::sd(y)
  }

  lstsq <- function(A) {
    # Perform the least squares fit
    #beta <- try(solve(t(A) %*% A) %*% t(A) %*% y, silent = FALSE)
    tryCatch({
      # Least squares solver
      beta <- solve(t(A) %*% A) %*% t(A) %*% y
      #utils::capture.output(bb <- fit(n_segments = 3))
    }, error = function(e) {
      ssr <<- Inf
      return(ssr)
    })
    y_hat <- A %*% beta
    e <- y_hat - y
    ssr <<- sum(e^2)
    return(ssr)
  }


  fit <- function(n_segments) {
    # Define the function to minimize
    min_function <- fit_with_breaks_opt

    # Store the number of line segments and number of parameters
    n_segments <- as.integer(n_segments)
    n_parameters <- n_segments + 1

    # Calculate the number of variables to solve for
    nVar <- n_segments - 1

    # Initiate the bounds of the optimization
    bounds <- matrix(0, nrow = nVar, ncol = 2)
    bounds[, 1] <- min(x)
    bounds[, 2] <- max(x)

    # Run the optimization
    #init <- lhs::randomLHS(100, nVar) * (max(x) - min(x)) + min(x)
    #init <- lapply(1:nrow(init), function(j) {sort(init[j,])}) %>% do.call("rbind", .)

    res <- DEoptim::DEoptim(
      min_function, lower = bounds[, 1], upper = bounds[, 2],
      control = DEoptim::DEoptim.control(
        VTR = 0,
        NP = 100,
        itermax = n_trials,
        reltol = 1e-3,
        CR = 0.7,
        strategy = 2,
        F = 0.8,
        steptol = 100
        #initialpop = init
      )
    )

    return(sort(res$optim$bestmem))
  }

  fit_with_breaks_opt <- function(breaks) {
    # Ensure necessary attributes are initialized
    if(!all(((biPOD:::bp_to_groups(dplyr::tibble(time=x, count=y), break_points = breaks) %>% table()) >= avg_points_per_window))) return(Inf)
    A <- assemble_regression_matrix(c(min(x), breaks, max(x)))

    # Try to solve the regression problem
    tryCatch({
      # Least squares solver
      ssr <<- lstsq(A)
      if (is.null(ssr)) {
        ssr <<- Inf
      }
    }, error = function(e) {
      # The computation could not converge
      ssr <<- Inf
    })

    return(ssr)
  }

  assemble_regression_matrix <- function(breaks) {
    # Ensure breaks is a numeric vector
    breaks <- as.numeric(breaks)
    # Sort the breaks and store them
    #breaks <- sort(breaks)
    fit_breaks <- breaks
    n_segments <- length(breaks) - 1

    # Assemble the regression matrix
    A_list <- list(rep(1, length(x)))

    A_list[[length(A_list) + 1]] <- x - fit_breaks[1]
    if ((n_segments - 1) >= 1) {
      for (i in 1:(n_segments - 1)) {
        A_list[[length(A_list) + 1]] <- ifelse(x > fit_breaks[i + 1], x - fit_breaks[i + 1], 0)
      }
    }

    A <- do.call(cbind, A_list)
    return(A)
  }

  max_breakpoints = min(max_breakpoints, as.integer(length(x) / avg_points_per_window))
  available_breakpoints <- 1:max_breakpoints
  message("Intializing breakpoints")
  proposed_breakpoints <- parallel::mclapply(available_breakpoints, FUN = function(n_breakpoints) {
    print(n_breakpoints)
    convergence <<- FALSE
    iter <<- 1
    bb <<- NULL

    tmp <- utils::capture.output(
      while (!convergence & iter < n_trials) {
        tryCatch({
          # Least squares solver
          bb <<- fit(n_segments = n_breakpoints + 1)
          convergence <<- TRUE
        }, error = function(e) {
          # The computation could not converge
          print("error")
          convergence <<- FALSE
          iter <<- iter + 1
        })
      }
    )

    dplyr::tibble(n_breakpoints = n_breakpoints, best_bp = list(bb), convergence = convergence)
  }, mc.cores = parallel::detectCores()) %>% do.call("bind_rows", .) %>% dplyr::distinct()

  for (j in 1:nrow(proposed_breakpoints)) {
    if (!all((biPOD:::bp_to_groups(dplyr::tibble(time=x, count=y), unlist(proposed_breakpoints[j,]$best_bp)) %>% table()) >= avg_points_per_window)) {
      proposed_breakpoints[j,]$convergence <- F
    }
  }

  proposed_breakpoints <- proposed_breakpoints %>% dplyr::filter(convergence)

  if (nrow(proposed_breakpoints) == 0) {
    message("Zero models with breakpoints has been found")
    return(list(best_bp=NULL, best_fit=NULL))
  }

  m <- biPOD:::get_model("fit_breakpoints")
  message("Breakpoints optimization")
  fits <- list()
  plots <- list()
  proposed_breakpoints$idx <- c(1:nrow(proposed_breakpoints)) + 1


  loos <- lapply(0:nrow(proposed_breakpoints), function(j) {
    #print(j)
    if (j == 0) {
      bp = array(0, dim = c(0))
    } else {
      bp <- sort(unlist(proposed_breakpoints[j,2]))
    }

    input_data <- list(
      S = length(x),
      G = length(bp),
      N = y,
      T = x,
      b_prior = bp,
      b = bp,
      #sigma_changepoints = (max(x) - min(x)) / length(bp)
      sigma_changepoints = (max(x) - min(x))
    )

    tmp <- utils::capture.output(
      suppressMessages(
        f <- m$sample(input_data, parallel_chains = 4)
      )
    )

    fits[[j+1]] <<- f
    suppressWarnings(loo <- f$loo())

    repetitions <- f$draws("y_rep", format = "matrix")
    means <- lapply(1:ncol(repetitions), function(j) {mean(repetitions[,j])}) %>% unlist()
    sds <- lapply(1:ncol(repetitions), function(j) {sd(repetitions[,j])}) %>% unlist()

    # plots[[j+1]] <<- ggplot2::ggplot() +
    #   ggplot2::geom_pointrange(dplyr::tibble(x=x, means=means, sds = sds), mapping = ggplot2::aes(x=.data$x, y=.data$means, ymin=.data$means-.data$sds, ymax=.data$means+.data$sds)) +
    #   ggplot2::geom_point(dplyr::tibble(x=x, y=y), mapping=ggplot2::aes(x=.data$x, y=.data$y), col="red") +
    #   ggplot2::geom_vline(xintercept = bp) +
    #   ggplot2::ggtitle(max(f$lp()))


    k = 1 + (1 + length(bp)) #+ length(bp)
    n = length(x)
    bic = k * log(n) - 2 * max(f$lp())
    #print(bic)

    return(bic)

    #return(loo)

    # k = 1 + (1 + length(bp)) #+ length(bp)
    # n = length(x)
    # bic = k * log(n) - 2 * max(f$lp())
    # #median(f$lp())
    # bic
  }) %>% unlist()

  #proposed_breakpoints
  #ggpubr::ggarrange(plotlist = plots)

  if (length(loos) == 1) {
    message("Zero models with breakpoints has been found")
    return(list(best_bp=NULL, best_fit=NULL))
  }

  message("Choosing optimal breakpoints")

  bic_comp <- dplyr::tibble(bic = loos, n_breakpoints = 1:length(loos) - 1)
  best_j <- bic_comp %>%
    dplyr::filter(bic == min(bic)) %>%
    dplyr::pull(n_breakpoints) %>%
    min()

  # suppressWarnings(loo_comp <- loo::loo_compare(loos))
  # loo_comp[,1] <- round(loo_comp[,1],1)
  #
  # loo_comp <- loo_comp %>% as_tibble() %>%
  #   dplyr::mutate(model = rownames(loo_comp)) %>%
  #   dplyr::mutate(j = as.numeric(stringr::str_replace(rownames(loo_comp), pattern = "model", replacement = "")) - 1) %>%
  #   dplyr::mutate(convergence = TRUE)
  #
  # loo_comp$convergence <- lapply(loo_comp$j, function(j) {
  #   if (j == 0) return(TRUE)
  #   all(loos[[j]]$diagnostics$pareto_k <= 1)
  # }) %>% unlist()
  #
  # loo_comp <- loo_comp %>% dplyr::filter(convergence) %>% dplyr::filter(elpd_diff == max(elpd_diff))
  # best_j <- min(loo_comp$j)

  best_fit <- fits[[best_j + 1]]
  best_fit <- biPOD:::convert_mcmc_fit_to_biPOD(best_fit)

  best_bp <- proposed_breakpoints %>%
    dplyr::filter(idx == best_j + 1) %>%
    dplyr::pull(best_bp) %>% unlist()

  if (norm) {
    x <- d$time
    best_bp <- best_bp * stats::sd(x) + mean(x)
  }

  return(list(best_bp=best_bp, best_fit=best_fit, plots=plots))
}
