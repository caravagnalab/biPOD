#' Print for class \code{'bipod'}.
#'
#' @param x An obj of class \code{'bipod'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @exportS3Method print bipod
#' @export print.bipod
print.bipod = function(x, ...) {
  stopifnot(inherits(x, "bipod"))

  nobs <- nrow(x$counts)
  nwindows <- length(unique(x$counts$group))

  cli::cli_rule(
    paste(
      crayon::bgYellow(
        crayon::black(
          paste0("[ biPOD ]", x$sample))),
      '{.field {nobs}} observations divided in {.field {nwindows}} time windows.'
    )
  )

  task1 = !is.null(x$fit)
  task2 = !is.null(x$breakpoints_fit)
  task3 = !is.null(x$two_pop_fit)

  # task 2 ####
  if (task2) {
    cat('\n')

    if (!is.null(x$metadata$breakpoints)) {

      rhats <- x$breakpoints_fit$rhat %>% unlist()
      mean_rhat <- c(rhats[grepl("b[", x$breakpoints_fit$parameters, fixed = TRUE)] %>% mean() %>% round(2))

      if (mean_rhat <= 1.1) {
        pstatus = function()
          "{crayon::bgGreen(crayon::black(\" PASS \"))}"
      } else {
        pstatus = function()
          "{crayon::bgRed(crayon::black(\" FAIL \"))}"
      }

      cli::cli_h1(
        paste0(
          "Break-points inference ",
          pstatus(),
          " Mean rhat = {crayon::green(mean_rhat)}."
        )
      )

      cli::cli_alert_info(paste0(" Number of breakpoints inferred : ", crayon::blue(length(x$metadata$breakpoints)), ""))

      pars <- lapply(c("b["), function(s) {
        x$breakpoints_fit$parameters[grepl(s, x$breakpoints_fit$parameters, fixed = TRUE)]
      }) %>% unlist() %>% unique()

      par_tibble <- lapply(pars, function(p) {
        draws <- x$breakpoints_fit$draws[,grepl(p, colnames(x$breakpoints_fit$draws), fixed = T)] %>% as.vector() %>% unlist()
        dplyr::tibble(
          Parameter = p,
          Mean = mean(draws),
          Sd = stats::sd(draws),
          p05 = stats::quantile(draws, .05),
          p50 = stats::quantile(draws, .50),
          p95 = stats::quantile(draws, .95),
          Rhat = x$breakpoints_fit$rhat[p] %>% as.numeric())
      }) %>% do.call(dplyr::bind_rows, .)
      print(par_tibble)
    } else {
      cli::cli_h1(
        paste0(
          "Task 0 ",
          " Break-points inference"
        )
      )

      cli::cli_alert_info(paste0(" Number of breakpoints inferred : ", crayon::blue(length(x$metadata$breakpoints)), ""))
    }
  }

  # task 3 ####
  if (task3) {
    cat('\n')

    if (x$metadata$status == "PASS") {
      pstatus = function()
        "{crayon::bgGreen(crayon::black(\" PASS \"))}"
    } else {
      pstatus = function()
        "{crayon::bgRed(crayon::black(\" FAIL \"))}"
    }

    mean_rhat <- x$two_pop_fit$rhat %>% unlist() %>% mean() %>% round(2)

    cli::cli_h1(
      paste0(
        "Two population inference ",
        pstatus(),
        " Using {crayon::bold(x$metadata$sampling)}. Mean rhat = {crayon::green(mean_rhat)}."
      )
    )

    x$two_pop_fit

    if (x$metadata$factor_size != 1) {
      cli::cli_alert_info(paste0(" Scale factor : ", crayon::blue(x$metadata$factor_size), ". Instant of birth might be incorrect."))
    } else {
      cli::cli_alert_info(paste0(" Scale factor : ", crayon::blue(x$metadata$factor_size), ""))
    }
    cat("\n")
    cli::cli_alert_info(" Inferred parameters")

    par_tibble <- lapply(x$two_pop_fit$parameters, function(p) {
      draws <- x$two_pop_fit$draws[,grepl(p, colnames(x$two_pop_fit$draws), fixed = T)] %>% as.vector() %>% unlist()
      dplyr::tibble(Parameter = p, Mean = mean(draws), Sd = stats::sd(draws), p05 = stats::quantile(draws, .05), p50 = stats::quantile(draws, .50), p95 = stats::quantile(draws, .95))
    }) %>% do.call(dplyr::bind_rows, .)
    print(par_tibble)
  }


  # task 1 ####
  if (task1) {
    cat('\n')

    if (x$metadata$status == "PASS") {
      pstatus = function()
        "{crayon::bgGreen(crayon::black(\" PASS \"))}"
    } else {
      pstatus = function()
        "{crayon::bgRed(crayon::black(\" FAIL \"))}"
    }

    mean_rhat <- x$fit$rhat %>% unlist() %>% mean() %>% round(2)

    cli::cli_h1(
      paste0(
        "Single population inference ",
        pstatus(),
        " Using {crayon::bold(x$metadata$sampling)}. Mean rhat = {crayon::green(mean_rhat)}."
      )
    )

    cli::cli_alert_info(paste0(" Growth pattern : ", crayon::blue(x$metadata$growth_type), ""))
    if (x$metadata$t0_inferred) {
      cli::cli_alert_info(paste0(" Instant of birth : ", crayon::blue("inferred"), ""))
    } else {
      cli::cli_alert_info(paste0(" Instant of birth : ", crayon::blue("not inferred"), ""))
    }
    if (x$metadata$factor_size != 1) {
      cli::cli_alert_info(paste0(" Scale factor : ", crayon::blue(x$metadata$factor_size), ". Instant of birth might be incorrect."))
    } else {
      cli::cli_alert_info(paste0(" Scale factor : ", crayon::blue(x$metadata$factor_size), ""))
    }
    cat("\n")
    cli::cli_alert_info(" Inferred parameters")

    par_tibble <- lapply(x$fit$parameters, function(p) {
      draws <- x$fit$draws[,grepl(p, colnames(x$fit$draws), fixed = T)] %>% as.vector() %>% unlist()
      dplyr::tibble(Parameter = p, Mean = mean(draws), Sd = stats::sd(draws), p05 = stats::quantile(draws, .05), p50 = stats::quantile(draws, .50), p95 = stats::quantile(draws, .95))
    }) %>% do.call(dplyr::bind_rows, .)
    print(par_tibble)
  }
}
