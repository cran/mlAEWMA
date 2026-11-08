# Generic engine shared by mlaewma_neural(), mlaewma_svm(), mlaewma_tree() and
# mlaewma_dqr().

.prepare_subgroups <- function(data, subgroup_size) {
  if (is.null(data) || length(data) == 0) {
    stop("`data` must not be empty.", call. = FALSE)
  }
  if (is.data.frame(data)) {
    if (!all(vapply(data, is.numeric, logical(1)))) {
      stop("`data` must be numeric.", call. = FALSE)
    }
    data <- as.matrix(data)
  }
  if (is.matrix(data)) {
    if (!is.numeric(data)) stop("`data` must be numeric.", call. = FALSE)
    if (ncol(data) != 5) {
      stop("`data` has ", ncol(data), " column(s) but subgroups of size 5 ",
           "are required: give one subgroup per row (5 columns), or a vector.",
           call. = FALSE)
    }
    application_data <- data
  } else if ((is.numeric(data) || is.logical(data)) && is.null(dim(data))) {
    # plain vector (also a named vector or a time series)
    data <- as.numeric(data)
    n_complete <- floor(length(data) / subgroup_size)
    if (n_complete < 1) {
      stop("Not enough observations to form one complete subgroup.",
           call. = FALSE)
    }
    # an incomplete final subgroup is ignored
    data <- data[seq_len(n_complete * subgroup_size)]
    application_data <- matrix(data, ncol = subgroup_size, byrow = TRUE)
  } else {
    stop("`data` must be a numeric vector, matrix or data frame.",
         call. = FALSE)
  }
  if (any(!is.finite(application_data))) {
    stop("`data` must contain no missing or infinite values.", call. = FALSE)
  }
  application_data
}

.check_h <- function(h, name = "h") {
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0) {
    stop("`", name, "` must be a single finite number > 0.", call. = FALSE)
  }
  invisible(h)
}

# Control limits for method = "all": NULL -> defaults, named -> matched
# case-insensitively, unnamed -> positional.
.resolve_h_all <- function(h, methods) {
  if (is.null(h)) {
    return(as.list(.aewma_default_h[methods]))
  }
  if (is.list(h)) h <- unlist(h)
  if (!is.numeric(h)) stop("`h` must be numeric.", call. = FALSE)
  for (i in seq_along(h)) .check_h(h[[i]], "h")
  if (!is.null(names(h))) {
    lower <- tolower(names(h))
    if (!setequal(lower, methods)) {
      stop("For method = 'all', a named h vector must have names ",
           paste(sQuote(methods, FALSE), collapse = ", "),
           " (case-insensitive). Got: ", paste(names(h), collapse = ", "),
           call. = FALSE)
    }
    names(h) <- lower
    return(as.list(h[methods]))
  }
  if (length(h) != length(methods)) {
    stop("For method = 'all', unnamed h must contain exactly ",
         length(methods), " values.", call. = FALSE)
  }
  stats::setNames(as.list(h), methods)
}

.find_h <- function(method, target_arl0, arl_tolerance, h_search_nsim,
                    h_validation_nsim, max_rl, seed, n_cores, n_streams) {
  lower <- 0.01
  upper <- .aewma_h_upper[[method]]

  for (iteration in 1:20) {
    middle <- (lower + upper) / 2
    result <- .arl(method, middle, 1, h_search_nsim, max_rl, seed, n_cores,
                   n_streams)
    message(toupper(method), " iteration = ", iteration,
            " h = ", round(middle, 4),
            " ARL0 = ", round(result["ARL"], 2))

    if (abs(result["ARL"] - target_arl0) <= arl_tolerance) break
    if (result["ARL"] < target_arl0) lower <- middle else upper <- middle
  }

  final_validation <- .arl(method, middle, 1, h_validation_nsim, max_rl,
                           seed + 999, n_cores, n_streams)
  c(h = middle,
    ARL = unname(final_validation["ARL"]),
    SDRL = unname(final_validation["SDRL"]),
    ARL_SE = unname(final_validation["ARL_SE"]),
    Censor_Rate = unname(final_validation["Censor_Rate"]))
}

.perf_row <- function(method, shift, perf) {
  data.frame(
    Method = toupper(method),
    Shift = shift,
    ARL = perf["ARL"],
    SDRL = perf["SDRL"],
    ARL_SE = perf["ARL_SE"],
    Censor_Rate = perf["Censor_Rate"]
  )
}

# Best method per shift (all-methods output).
.best_per_shift <- function(performance, shifts, target_arl0) {
  best <- data.frame()
  for (shift in shifts) {
    temp <- performance[performance$Shift == shift, ]

    if (any(temp$Censor_Rate > 0.01)) {
      best_row <- NA_integer_
      status <- "Censored - increase max_rl"
    } else if (shift == 1) {
      best_row <- which.min(abs(temp$ARL - target_arl0))
      status <- "Closest to target ARL0"
    } else {
      best_row <- which.min(temp$ARL)
      status <- "Fastest detection by ARL"
    }

    if (is.na(best_row)) {
      best_method <- NA_character_
      best_arl <- best_sdr <- best_arl_se <- NA_real_
    } else {
      best_method <- temp$Method[best_row]
      best_arl <- temp$ARL[best_row]
      best_sdr <- temp$SDRL[best_row]
      best_arl_se <- temp$ARL_SE[best_row]
    }

    best <- rbind(best, data.frame(
      Shift = shift, Best_Method = best_method, Best_ARL = best_arl,
      Best_SDRL = best_sdr, Best_ARL_SE = best_arl_se,
      Comparison_Status = status
    ))
  }
  best
}

# Average-rank summary over shifts where all methods are uncensored.
.rank_summary <- function(performance, method_names) {
  rank_data <- performance[performance$Shift != 1, ]
  if (nrow(rank_data) == 0) return(NULL)

  valid_shifts <- numeric(0)
  for (shift_value in unique(rank_data$Shift)) {
    temp_shift <- rank_data[rank_data$Shift == shift_value, ]
    if (nrow(temp_shift) == length(method_names) &&
          all(method_names %in% temp_shift$Method) &&
          all(temp_shift$Censor_Rate <= 0.01)) {
      valid_shifts <- c(valid_shifts, shift_value)
    }
  }

  if (length(valid_shifts) == 0) {
    return(data.frame(Method = character(0), Mean_Rank = numeric(0),
                      stringsAsFactors = FALSE))
  }

  rank_rows <- data.frame(Method = character(0), Shift = numeric(0),
                          Rank = numeric(0), stringsAsFactors = FALSE)
  for (shift_value in valid_shifts) {
    temp_shift <- rank_data[rank_data$Shift == shift_value, ]
    temp_shift <- temp_shift[match(method_names, temp_shift$Method), ]
    rank_rows <- rbind(rank_rows, data.frame(
      Method = temp_shift$Method, Shift = shift_value,
      Rank = rank(temp_shift$ARL, ties.method = "average"),
      stringsAsFactors = FALSE
    ))
  }

  rank_summary <- stats::aggregate(Rank ~ Method, data = rank_rows, FUN = mean)
  names(rank_summary)[names(rank_summary) == "Rank"] <- "Mean_Rank"
  rank_summary <- rank_summary[match(method_names, rank_summary$Method), ]
  rownames(rank_summary) <- NULL
  rank_summary
}

.aewma_family <- function(family, data, method, h, sigma0, subgroup_size,
                          find_h, target_arl0, h_search_nsim,
                          h_validation_nsim, arl_tolerance, shifts, shift_nsim,
                          shift_max_rl, nsim, max_rl, seed, n_cores,
                          n_streams) {
  methods <- .aewma_families[[family]]
  lagged <- family == "dqr"
  fun <- paste0("mlaewma_", family, "()")
  # `all` is the same as `dqr` for the one-model DQR family
  if (lagged && is.character(method) && identical(tolower(method), "all")) {
    method <- "dqr"
  }
  method <- .match_method(method, c(methods, if (!lagged) "all"), fun)

  .check_subgroup_size(subgroup_size)
  if (!is.null(sigma0) && (!is.numeric(sigma0) || length(sigma0) != 1L ||
                             !is.finite(sigma0) || sigma0 <= 0)) {
    stop("`sigma0` must be a single positive number.", call. = FALSE)
  }
  .check_count(n_cores, "n_cores")
  .check_count(n_streams, "n_streams")
  .check_seed(seed)
  if (!(is.logical(find_h) && length(find_h) == 1L && !is.na(find_h))) {
    stop("`find_h` must be TRUE or FALSE.", call. = FALSE)
  }
  if (find_h) {
    .check_pos(target_arl0, "target_arl0")
    .check_pos(arl_tolerance, "arl_tolerance")
  }
  if (!is.null(shifts)) .check_shifts(shifts)

  application_data <- .prepare_subgroups(data, subgroup_size)
  if (is.null(sigma0)) {
    sigma0 <- stats::sd(as.numeric(application_data))
    if (!is.finite(sigma0) || sigma0 <= 0) {
      stop("`sigma0` could not be estimated because the data have no ",
           "variation; supply `sigma0`.", call. = FALSE)
    }
  }

  # Castagliola (2005) transformation of the subgroup variance.
  a <- .AEWMA_A - 2 * .AEWMA_B * log(sigma0)
  c_value <- .AEWMA_C * sigma0^2
  Mx <- apply(application_data, 1, stats::var)
  D <- a + .AEWMA_B * log(Mx + c_value)
  if (any(!is.finite(D))) {
    stop("The transformed statistic is not finite; check the scale of ",
         "`data` and `sigma0`.", call. = FALSE)
  }

  theta_fun <- if (lagged) .theta_fun_dqr() else NULL

  run_method <- function(m, user_h = NULL) {
    calibration <- NULL

    if (isTRUE(find_h)) {
      calibration <- .find_h(m, target_arl0, arl_tolerance, h_search_nsim,
                             h_validation_nsim, max_rl, seed, n_cores,
                             n_streams)
      h_used <- as.numeric(calibration["h"])
    } else if (!is.null(user_h)) {
      .check_h(user_h, "h")
      h_used <- as.numeric(user_h)
    } else {
      h_used <- .aewma_default_h[[m]]
    }

    if (lagged) {
      theta <- numeric(length(D))
      D_lag1 <- 0
      D_lag2 <- 0
      for (i in seq_along(D)) {
        if (i == 1) {
          xt_lag1 <- 0
          xt_lag2 <- 0
        } else if (i == 2) {
          xt_lag1 <- D_lag1
          xt_lag2 <- 0
        } else {
          xt_lag1 <- D_lag1
          xt_lag2 <- D_lag2
        }
        theta[i] <- theta_fun(D[i], xt_lag1, xt_lag2)
        D_lag2 <- D_lag1
        D_lag1 <- D[i]
      }
    } else {
      theta <- .theta_fun(m)(D)
    }

    AEWMA <- numeric(length(D))
    AEWMA[1] <- theta[1] * D[1]
    if (length(D) > 1) {
      for (i in 2:length(D)) {
        AEWMA[i] <- theta[i] * D[i] + (1 - theta[i]) * AEWMA[i - 1]
      }
    }

    signal_points <- which(abs(AEWMA) > h_used)
    first_signal <- if (length(signal_points) == 0) NA_integer_
                    else signal_points[1]

    list(
      method = toupper(m),
      h = h_used,
      first_signal = first_signal,
      results = data.frame(
        Subgroup = seq_along(D), Variance = Mx, D = D, Theta = theta,
        AEWMA = AEWMA, h = h_used, Exceeds_Limit = abs(AEWMA) > h_used
      ),
      calibration = calibration
    )
  }

  shift_arl <- function(m, h_value, shift, shift_index) {
    .arl(m, h_value, shift, shift_nsim, shift_max_rl,
         seed + 1000L * shift_index, n_cores, n_streams)
  }

  # ---- one method -------------------------------------------------------
  if (method != "all") {
    result <- run_method(method, user_h = h)
    if (!is.null(shifts)) {
      performance <- data.frame()
      for (shift_index in seq_along(shifts)) {
        shift <- shifts[shift_index]
        temp <- shift_arl(method, result$h, shift, shift_index)
        performance <- rbind(performance, .perf_row(method, shift, temp))
      }
      result$performance <- performance
    }
    return(result)
  }

  # ---- all methods of the family -----------------------------------------
  h_list <- .resolve_h_all(h, methods)
  results <- lapply(methods, function(m) run_method(m, user_h = h_list[[m]]))
  names(results) <- toupper(methods)

  comparison <- data.frame(
    Method = toupper(methods),
    h = vapply(results, function(r) r$h, numeric(1), USE.NAMES = FALSE),
    First_Signal = vapply(results, function(r) r$first_signal, integer(1),
                          USE.NAMES = FALSE)
  )

  performance <- NULL
  best <- NULL
  rank_summary <- NULL
  calibration_results <- lapply(results, function(r) r$calibration)

  if (!is.null(shifts)) {
    performance <- data.frame()
    for (shift_index in seq_along(shifts)) {
      shift <- shifts[shift_index]
      rows <- lapply(methods, function(m) {
        .perf_row(m, shift,
                  shift_arl(m, results[[toupper(m)]]$h, shift, shift_index))
      })
      performance <- do.call(rbind, c(list(performance), rows))
    }
    best <- .best_per_shift(performance, shifts, target_arl0)
    rank_summary <- .rank_summary(performance, toupper(methods))
  }

  c(results, list(comparison = comparison, performance = performance,
                  best = best, rank_summary = rank_summary,
                  calibration = calibration_results))
}
