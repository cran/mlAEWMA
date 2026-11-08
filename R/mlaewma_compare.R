#' Compare all nine AEWMA charts
#'
#' Runs the neural network (`nn`, `qrnn`, `mcqrnn`), SVM (`svr`, `svqr`,
#' `rsvqr`), tree (`dt`, `gb`) and DQR (`dqr`) charts on the same data and
#' collects them in one table. With `shifts`, the out-of-control average run
#' lengths are simulated and the models are ranked.
#'
#' @inheritParams mlaewma_family
#' @param h `NULL` for the calibrated defaults, otherwise a named vector or
#'   list that contains all nine model names (`nn`, `qrnn`, `mcqrnn`, `svr`,
#'   `svqr`, `rsvqr`, `dt`, `gb`, `dqr`; case-insensitive).
#' @param n_cores Number of worker processes for simulation; `NULL` means 1.
#' @param source_dir Not used any more; kept so that older calls still work.
#' @return A list with `comparison` (one row per model: `Family`, `Method`,
#'   `h`, `First_Signal`), `models` (the nine individual results),
#'   `performance`, `rank_summary` and `family_results` (the four family
#'   results).
#' @examples
#' data(patient)
#' mlaewma_compare(patient, sigma0 = 4.05)$comparison
#' @export
mlaewma_compare <- function(data, sigma0 = NULL, subgroup_size = 5, h = NULL,
                          find_h = FALSE, target_arl0 = 370,
                          h_search_nsim = 3000, h_validation_nsim = 20000,
                          arl_tolerance = 5, shifts = NULL,
                          shift_nsim = 10000, shift_max_rl = 10000,
                          nsim = 500, max_rl = 10000, seed = 42,
                          n_cores = NULL, source_dir = ".",
                          n_streams = 3L) {
  if (!identical(source_dir, ".")) {
    warning("`source_dir` is not used any more.", call. = FALSE)
  }
  if (is.null(n_cores)) n_cores <- 1L

  if (!is.null(h)) {
    if (is.null(names(h))) {
      stop("h must be a named vector/list; use names(): ",
           paste(.ALL_METHODS, collapse = ", "), call. = FALSE)
    }
    h_names_lower <- tolower(names(h))
    if (!setequal(h_names_lower, .ALL_METHODS)) {
      stop("h must contain exactly these 9 names (case-insensitive): ",
           paste(.ALL_METHODS, collapse = ", "), ". Got: ",
           paste(names(h), collapse = ", "), call. = FALSE)
    }
    names(h) <- h_names_lower
  }
  h_for <- function(method_names) if (is.null(h)) NULL else h[method_names]

  args <- list(
    sigma0 = sigma0, subgroup_size = subgroup_size, find_h = find_h,
    target_arl0 = target_arl0, h_search_nsim = h_search_nsim,
    h_validation_nsim = h_validation_nsim, arl_tolerance = arl_tolerance,
    shifts = shifts, shift_nsim = shift_nsim, shift_max_rl = shift_max_rl,
    nsim = nsim, max_rl = max_rl, seed = seed, n_cores = n_cores,
    n_streams = n_streams
  )

  neural_result <- do.call(mlaewma_neural, c(
    list(data = data, method = "all", h = h_for(c("nn", "qrnn", "mcqrnn"))),
    args))
  svm_result <- do.call(mlaewma_svm, c(
    list(data = data, method = "all", h = h_for(c("svr", "svqr", "rsvqr"))),
    args))
  tree_result <- do.call(mlaewma_tree, c(
    list(data = data, method = "all", h = h_for(c("dt", "gb"))), args))
  dqr_h <- h_for("dqr")
  dqr_result <- do.call(mlaewma_dqr, c(
    list(data = data, method = "dqr",
         h = if (is.null(dqr_h)) NULL else as.numeric(dqr_h)), args))

  models <- list(
    NN = neural_result$NN, QRNN = neural_result$QRNN,
    MCQRNN = neural_result$MCQRNN,
    SVR = svm_result$SVR, SVQR = svm_result$SVQR, RSVQR = svm_result$RSVQR,
    DT = tree_result$DT, GB = tree_result$GB,
    DQR = dqr_result
  )
  method_display_names <- c("NN", "QRNN", "MCQRNN", "SVR", "SVQR", "RSVQR",
                            "DT", "GB", "DQR")
  family_of <- c(
    NN = "Neural", QRNN = "Neural", MCQRNN = "Neural",
    SVR = "SVM", SVQR = "SVM", RSVQR = "SVM",
    DT = "Tree", GB = "Tree", DQR = "DQR"
  )

  comparison <- data.frame(
    Family = unname(family_of[method_display_names]),
    Method = method_display_names,
    h = vapply(models, function(m) as.numeric(m$h), numeric(1)),
    First_Signal = vapply(
      models,
      function(m) if (is.na(m$first_signal)) NA_integer_ else m$first_signal,
      integer(1)
    ),
    row.names = NULL
  )

  performance <- NULL
  rank_summary <- NULL
  if (!is.null(shifts)) {
    performance <- rbind(neural_result$performance, svm_result$performance,
                         tree_result$performance, dqr_result$performance)
    rownames(performance) <- NULL

    if (nrow(performance) > 0) {
      # the shift 1 rows measure ARL0 and are not ranked
      rank_summary <- .rank_summary(performance, method_display_names)
      if (!is.null(rank_summary) && nrow(rank_summary) > 0) {
        rank_summary <- rank_summary[order(rank_summary$Mean_Rank), ]
        rownames(rank_summary) <- NULL
      }
    }
  }

  list(
    comparison = comparison,
    models = models,
    performance = performance,
    rank_summary = rank_summary,
    family_results = list(neural = neural_result, svm = svm_result,
                          tree = tree_result, dqr = dqr_result)
  )
}
