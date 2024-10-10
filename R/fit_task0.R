#' Fit breakpoints to bipod object
#'
#' @param x A bipod object.
#' @param norm Logical value indicating whether to normalize the data.
#'  If TRUE, the time and count data are standardized before fitting the model. (default is TRUE)
#' @param n_trials Integer specifying the number of trials for the optimization algorithm.
#'  This controls the number of iterations used to fit the breakpoints. (default is 500)
#' @param avg_points_per_window Integer specifying the average number of data points per segment.
#'  This parameter influences the granularity of the segments when searching for breakpoints. (default is 3)
#' @param available_changepoints Integer vector specifying the range of available changepoints.
#'  These values represent the possible number of breakpoints to be considered during model fitting. (default is 0:5)
#' @param model_selection Character string specifying the model selection criterion.
#'  Options include 'LOO' (Leave-One-Out cross-validation), 'AIC' (Akaike Information Criterion), and 'BIC' (Bayesian Information Criterion).  (default is "LOO")
#' @param n_core Integer specifying the number of CPU cores to use for parallel processing. (default is 4)
#'
#' @return The input bipod object with an added 'breakpoints_fit' slot containing the fitted model for the breakpoints.
#' @export
fit_breakpoints = function(
    x,
    norm=TRUE,
    n_trials=500,
    avg_points_per_window = 3,
    available_changepoints = c(0:5),
    model_selection = "LOO",
    n_core = 4
) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  #if (!(factor_size > 0)) stop("factor_size must be positive")
  d <- x$counts

  res <- find_breakpoints(
    d,
    avg_points_per_window=avg_points_per_window,
    available_changepoints=available_changepoints,
    norm=norm,
    n_trials=n_trials,
    model_selection=model_selection,
    n_core=n_core
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
  x$breakpoints_fit$normalization_pars <- res$normalization_pars

  if (!(is.null(best_bp))) {
    x$counts$group <- bp_to_groups(x$counts, x$metadata$breakpoints)
  }
  cli::cli_alert_info("Median of the inferred breakpoints have been succesfully stored.")

  x
}


find_breakpoints = function(d, avg_points_per_window, max_breakpoints, norm, n_trials,available_changepoints, model_selection, n_core) {
  x <- d$time
  y <- log(d$count)

  if (norm) {
    x <- (x - mean(x)) / stats::sd(x)
    y <- (y - mean(y)) / stats::sd(y)
  }

  lstsq <- function(A) {
    # Perform the least squares fit
    tryCatch({
      # Least squares solver
      beta <- solve(t(A) %*% A) %*% t(A) %*% y
    }, error = function(e) {
      ssr <<- Inf
      #ssr <- Inf
      return(ssr)
    })
    y_hat <- A %*% beta
    e <- y_hat - y
    ssr <- sum(e^2)
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
    if(!all(((bp_to_groups(dplyr::tibble(time=x, count=y), break_points = breaks) %>% table()) >= avg_points_per_window))) return(Inf)
    A <- assemble_regression_matrix(c(min(x), breaks, max(x)))

    # Try to solve the regression problem
    tryCatch({
      # Least squares solver
      ssr <- lstsq(A)
      if (is.null(ssr)) {
        return(Inf)
        #return(ssr)
      }
    }, error = function(e) {
      # The computation could not converge
      #ssr <<- Inf
      return(Inf)
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


  max_breakpoints = min(max(available_changepoints), as.integer(length(x) / avg_points_per_window))
  available_changepoints <- available_changepoints[available_changepoints <= max_breakpoints & available_changepoints > 0] %>% sort()
  if (!length(available_changepoints)) {
    cli::cli_alert_info("No breakpoints possible given the 'available_changepoints' and 'avg_points_per_window' passed as input!")
    cli::cli_alert_info("Model with zero breakpoints will be considered.")
    return(list(best_bp=NULL, best_fit=NULL))
  }

  cli::cli_alert_info("Intializing breakpoints...")
  proposed_breakpoints <- parallel::mclapply(available_changepoints, FUN = function(n_breakpoints) {
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
  }, mc.cores = min(parallel::detectCores(), n_core)) %>%
    do.call(dplyr::bind_rows, .) %>%
    dplyr::distinct()

  for (j in 1:nrow(proposed_breakpoints)) {
    if (!all((bp_to_groups(dplyr::tibble(time=x, count=y), unlist(proposed_breakpoints[j,]$best_bp)) %>% table()) >= avg_points_per_window)) {
      proposed_breakpoints[j,]$convergence <- F
    }
  }

  proposed_breakpoints <- proposed_breakpoints %>%
    dplyr::filter(.data$convergence == TRUE) %>%
    dplyr::mutate(idx = dplyr::row_number())

  if (nrow(proposed_breakpoints) == 0) {
    cli::cli_alert_info("Zero models with breakpoints has been found!")
    return(list(best_bp=NULL, best_fit=NULL))
  }

  m <- get_model("fit_breakpoints")
  cli::cli_alert_info("Breakpoints optimization...")
  fits <- list()
  #plots <- list()

  #proposed_breakpoints$idx <- c(1:nrow(proposed_breakpoints)) + 1

  criterion <- lapply(0:nrow(proposed_breakpoints), function(j) {
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
      sigma_changepoints = .01
      #sigma_changepoints = (max(x) - min(x)) / length(bp)
      #sigma_changepoints = (max(x) - min(x))
    )

    tmp <- utils::capture.output(
      suppressMessages(
        f <- m$sample(input_data, parallel_chains = 4)
      )
    )

    fits[[j+1]] <<- f

    if (model_selection == 'LOO') {
      suppressWarnings(loo <- f$loo())
      return(loo)
    } else if (model_selection == "BIC") {
      k = 2 + (1 + length(bp)) #+ length(bp)
      n = length(x)
      return(k * log(n) - 2 * max(f$lp()))
    }else if (model_selection == "AIC") {
      k = 2 + (1 + length(bp)) #+ length(bp)
      return(2 * k - 2 * max(f$lp()))
    } else {
      stop("model_selection parmater not recognised")
    }

    # repetitions <- f$draws("y_rep", format = "matrix")
    # means <- lapply(1:ncol(repetitions), function(j) {mean(repetitions[,j])}) %>% unlist()
    # sds <- lapply(1:ncol(repetitions), function(j) {sd(repetitions[,j])}) %>% unlist()

    # plots[[j+1]] <<- ggplot2::ggplot() +
    #   ggplot2::geom_pointrange(dplyr::tibble(x=x, means=means, sds = sds), mapping = ggplot2::aes(x=.data$x, y=.data$means, ymin=.data$means-.data$sds, ymax=.data$means+.data$sds)) +
    #   ggplot2::geom_point(dplyr::tibble(x=x, y=y), mapping=ggplot2::aes(x=.data$x, y=.data$y), col="red") +
    #   ggplot2::geom_vline(xintercept = bp) +
    #   ggplot2::ggtitle(max(f$lp()))


    #return(loo)
  })

  #proposed_breakpoints
  #ggpubr::ggarrange(plotlist = plots)

  if (length(criterion) == 1) {
    message("Zero models with breakpoints has been found")
    return(list(best_bp=NULL, best_fit=NULL))
  }

  cli::cli_alert_info("Choosing optimal breakpoints...")

  if (model_selection == 'LOO') {
    suppressWarnings(loo_comp <- loo::loo_compare(criterion))
    loo_comp[,1] <- round(loo_comp[,1],1)

    loo_comp <- loo_comp %>% dplyr::as_tibble() %>%
      dplyr::mutate(model = rownames(loo_comp)) %>%
      dplyr::mutate(j = as.numeric(stringr::str_replace(rownames(loo_comp), pattern = "model", replacement = ""))) %>%
      dplyr::mutate(convergence = TRUE)

    loo_comp$convergence <- lapply(loo_comp$j, function(j) {
      if (j == 0) return(TRUE)
      all(criterion[[j]]$diagnostics$pareto_k <= 1.1)
    }) %>% unlist()

    #loo_comp <- loo_comp %>% dplyr::filter(convergence) %>% dplyr::filter(elpd_diff == max(elpd_diff))
    loo_comp <- loo_comp %>% dplyr::filter(convergence) %>% dplyr::arrange(-as.numeric(.data$elpd_diff), -as.numeric(.data$se_diff))
    best_js <- loo_comp$j
  } else if (model_selection %in% c("AIC", "BIC")) {

    comp <- dplyr::tibble(value = unlist(criterion), n_breakpoints = c(0, proposed_breakpoints$n_breakpoints)) %>%
      dplyr::mutate(idx = dplyr::row_number())

    best_js <- comp %>%
      dplyr::arrange(.data$value) %>%
      dplyr::pull(.data$idx)

  } else {
    stop("model_selection parmater not recognised")
  }

  all_correct <- FALSE
  j_idx <- 1
  while (!all_correct) {
    best_j <- best_js[j_idx]

    best_bp <- proposed_breakpoints %>%
      dplyr::filter(.data$idx == (best_j - 1)) %>%
      dplyr::pull(.data$best_bp) %>% unlist()

    x <- d$time
    y <- log(d$count)

    if (norm) {
      x <- (x - mean(x)) / stats::sd(x)
      y <- (y - mean(y)) / stats::sd(y)
      #x <- d$time
      #best_bp <- best_bp * stats::sd(x) + mean(x)
    }

    # Extra fit
    m <- get_model("piecewise_changepoints")

    if (is.null(best_bp) | length(best_bp) == 0) {
      bp = array(0, dim = c(0))
    } else {
      bp <- sort(best_bp)
    }

    input_data <- list(
      S = length(x),
      G = length(bp),
      N = y,
      T = x,
      b_prior = bp,
      sigma_changepoints = .1
    )

    tmp <- utils::capture.output(
      suppressMessages(
        f <- m$sample(input_data, parallel_chains = 4)
      )
    )

    final_fit <- convert_mcmc_fit_to_biPOD(f, variational = F)

    if (norm) {
      #x <- d$time
      #final_fit$draws[grepl("b[", colnames(final_fit$draws), fixed = T)] <- final_fit$draws[grepl("b[", colnames(final_fit$draws), fixed = T)] * stats::sd(x) + mean(x)
      normalization_pars = list("mean" = mean(d$time), "sd" = stats::sd(d$time))
    } else {
      normalization_pars = NULL
    }

    if (length(bp)) {
      final_bp <- f$draws(variables = 'b', format = 'matrix') %>%
        dplyr::as_tibble() %>%
        dplyr::summarise_all(stats::median) %>%
        as.numeric()

      if (norm) {
        x <- d$time
        final_bp <- final_bp * stats::sd(x) + mean(x)
      }
    } else {
      final_bp <- NULL
    }

    # bring final_bp to true x values if very close to each other
    final_bp <- lapply(final_bp, function(f_bp) {
      if (any(abs(f_bp - x) <= (max(x) - min(x)) / 1000)) {
        return(x[abs(f_bp - x) <= (max(x) - min(x)) / 1000])
      } else {
        return(f_bp)
      }
    }) %>% unlist()

    all_correct <- all((bp_to_groups(dplyr::tibble(time=x, count=y), unlist(final_bp)) %>% table()) >= avg_points_per_window) & mean(unlist(final_fit$rhat)) <= 1.1
    j_idx <- j_idx + 1
  }

  return(list(best_bp=final_bp, best_fit=final_fit, normalization_pars=normalization_pars))
}
