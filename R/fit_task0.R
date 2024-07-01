#' Fit growth model to bipod object
#'
#' @param x a bipod object
#' @param norm .
#' @param n_trials .
#' @param min_points .
#' @param available_breakpoints .
#' @param constrain_bp_on_x .
#'
#' @return the input bipod object with an added 'breakpoints_fit' slot containing the fitted model for the breakpoints
#' @export
fit_breakpoints <- function(
    x,
    norm=F,
    n_trials=5000,
    min_points=3,
    available_breakpoints=c(1:5),
    constrain_bp_on_x=F
    ) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  #if (!(factor_size > 0)) stop("factor_size must be positive")
  d <- x$counts

  res <- find_breakpoints_v3(
    d,
    norm=norm,
    n_trials=n_trials,
    min_points=min_points,
    available_breakpoints=available_breakpoints,
    constrain_bp_on_x=constrain_bp_on_x
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

  # Write fit info
  # x$metadata$sampling <- sampling
  #x$metadata$factor_size <- factor_size
  # x$metadata$prior_K <- input_data$prior_K

  # Add median of breakpoints
  # n_changepoints <- length(input_data$changing_times_prior)
  # breakpoints_names <- lapply(1:n_changepoints, function(i) {
  #   paste0("changing_times[", i, "]")
  # }) %>% unlist()

  # if (best_res$J == 0) {
  #   median_breakpoints = NULL
  # } else {
  #   median_breakpoints <- best_fit$draws(variables = 'b', format = 'matrix') %>%
  #     dplyr::as_tibble() %>%
  #     dplyr::summarise_all(stats::median) %>%
  #     as.numeric()
  #
  #   median_breakpoints <- median_breakpoints + min(x$counts$time)
  # }

  x$metadata$breakpoints <- best_bp

  if (!(is.null(best_bp))) {
    x$counts$group <- bp_to_groups(x$counts, x$metadata$breakpoints)
  }

  if (!constrain_bp_on_x) {
    cli::cli_alert_success("Breakpoints have been inferred. Inspect the results using the {.field plot_breakpoints_posterior} function.")
  }
  cli::cli_alert_info("Median of the inferred breakpoints have been succesfully stored.")

  x
}

ind <- function(x, y) { return(as.numeric(x >= y)) }

# Function to calculate the expected mean
expected_mean <- function(x, q, s, b) {
  G <- length(s)
  res <- q + x * s[1]
  for (g in 2:G) {
    res <- res + (x - b[g-1]) * s[g] * ind(x, b[g-1])
  }
  return(res)
}

find_breakpoints_v3 <- function(d, norm=T, n_trials=1000, min_points=3, available_breakpoints=c(1:6), constrain_bp_on_x=F) {
  x <- d$time
  y <- log(d$count)

  if (norm) {
    x <- (x - mean(x)) / stats::sd(x)
    y <- (y - mean(y)) / stats::sd(y)
  }

  available_breakpoints <- available_breakpoints[available_breakpoints != 0]

  message("Initial proposals")
  proposed_breakpoints <- lapply(available_breakpoints, function(n_breakpoints) {
    min_rmse <- Inf
    best_starts <- NULL
    j_proposed <- 0
    j_iter <- 0
    while (j_proposed < n_trials & j_iter < n_breakpoints * n_trials) {
      j_iter <- j_iter + 1
      if (constrain_bp_on_x) {
        random_starts <- sample(x, n_breakpoints, replace = F)
      } else {
        random_starts <- lhs::randomLHS(1, n_breakpoints)
        random_starts <- random_starts * (max(x) - min(x)) + min(x)
      }

      bp <- sort(random_starts)
      n_per_window <- biPOD:::bp_to_groups(dplyr::tibble(time=x, count=y), bp) %>% table()
      if (any(n_per_window < min_points) | length(n_per_window) != ncol(random_starts) + 1) {next}

      # build design matrix
      n_params = n_breakpoints + 2
      X = matrix(0, nrow = length(x), ncol = n_params)
      X[,1] = 1
      X[,2] = x
      tmp <- lapply(1:ncol(random_starts), function(k) {
        X[,k+2] <<- ifelse(x > bp[k], x - bp[k], 0)
      })

      params <- c(solve(t(X) %*% X) %*% t(X) %*% y)
      ypred = expected_mean(x, params[1], params[2:length(params)], bp)
      rmse = sqrt(mean((y - ypred)**2))

      if (rmse < min_rmse) {
        min_rmse <- rmse
        best_starts <- bp
      }
      j_proposed <- j_proposed + 1
    }

    if (min_rmse < Inf) {
      return(dplyr::tibble(rmse = min_rmse, bp = list(best_starts), n_breakpoints=n_breakpoints))
    } else {
      return(NULL)
    }
  }) %>% do.call("bind_rows", .) %>% dplyr::distinct()

  if (nrow(proposed_breakpoints) == 0) {
    message("Zero models with breakpoints has been found")
    return(list(best_bp=NULL, best_fit=NULL))
  }

  #tmp <- utils::capture.output(suppressMessages(m <- cmdstanr::cmdstan_model("piecewise_fixed_breakpoints.stan")))

  if (constrain_bp_on_x == T) {
    m <- biPOD:::get_model("pw_lin_fixed_b")
  } else {
    m <- biPOD:::get_model("piecewise_changepoints")
  }

  message("Proposals' optimization")
  fits <- list()
  j = 0
  proposed_breakpoints$idx <- c(1:nrow(proposed_breakpoints)) + 1
  loos <- lapply(0:nrow(proposed_breakpoints), function(j) {
    if (j == 0) {
      bp = array(0, dim = c(0))
    } else {
      bp <- sort(unlist(proposed_breakpoints[j,2]))
    }

    if (constrain_bp_on_x == T) {
      input_data <- list(
        S = length(x),
        G = length(bp),
        N = y,
        T = x,
        b = bp
      )
    } else {
      input_data <- list(
        S = length(x),
        G = length(bp),
        N = y,
        T = x,
        b_prior = bp,
        sigma_changepoints = 1
      )
    }

    tmp <- utils::capture.output(
      suppressMessages(
        f <- m$sample(input_data, parallel_chains = 4)
      )
    )

    suppressWarnings(loo <- f$loo())
    fits[[j+1]] <<- f
    loo
  })

  if (length(loos) == 1) {
    message("Zero models with breakpoints has been found")
    return(list(best_bp=NULL, best_fit=NULL))
  }

  suppressWarnings(loo_comp <- loo::loo_compare(loos))
  best_j <- as.numeric(stringr::str_replace(rownames(loo_comp)[1], pattern = "model", replacement = ""))

  if (constrain_bp_on_x) {
    if (best_j == 1) { return(NULL) }
    best_bp <- proposed_breakpoints %>% dplyr::filter(idx == best_j) %>% pull(bp) %>% unlist() %>% sort()
    best_fit <- NULL
  } else {
    best_fit <- fits[[best_j]]
    if (best_j == 1) {
      best_bp = NULL
    } else {
      best_bp <- best_fit$draws(variables = 'b', format = 'matrix') %>%
        dplyr::as_tibble() %>%
        dplyr::summarise_all(stats::median) %>%
        as.numeric()
    }
    best_fit <- biPOD:::convert_mcmc_fit_to_biPOD(best_fit)
  }

  if (norm) {
    x <- d$time
    best_bp <- best_bp * stats::sd(x) + mean(x)
  }
  return(list(best_bp=best_bp, best_fit=best_fit))
}
