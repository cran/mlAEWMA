# Monte Carlo run-length simulation.
#
# Results depend on (seed, n_streams, nsim) only. The work is split into
# `n_streams` chunks; chunk i uses the i-th L'Ecuyer-CMRG stream derived from
# `seed`. `n_cores` only decides how many worker processes execute the
# chunks and never changes the numbers. n_streams = 1 uses a single
# Mersenne-Twister stream.

# One chart, memoryless models: subgroups are generated in blocks and theta is
# predicted for the whole block. If the chart signals inside a block, the RNG
# is rewound so the next chain continues from the same point as one-by-one
# generation would.
.chain_blocked <- function(predict_theta, h_value, shift, subgroup_size,
                           max_rl, block) {
  AEWMA <- 0
  done <- 0

  repeat {
    size <- min(block, max_rl - done)
    seed_before <- get(".Random.seed", envir = globalenv())

    x <- matrix(
      stats::rnorm(size * subgroup_size, mean = 0, sd = sqrt(shift)),
      ncol = subgroup_size, byrow = TRUE
    )
    Mx <- rowSums((x - rowMeans(x))^2) / (subgroup_size - 1)
    D <- .AEWMA_A + .AEWMA_B * log(Mx + .AEWMA_C)
    theta <- predict_theta(D)

    for (k in seq_len(size)) {
      AEWMA <- theta[k] * D[k] + (1 - theta[k]) * AEWMA
      if (abs(AEWMA) > h_value) {
        assign(".Random.seed", seed_before, envir = globalenv())
        invisible(stats::rnorm(k * subgroup_size, mean = 0, sd = sqrt(shift)))
        return(done + k)
      }
    }

    done <- done + size
    if (done >= max_rl) return(max_rl)
  }
}

# One chart, DQR: theta depends on the two previous D values, so subgroups
# are processed one at a time.
.chain_lagged <- function(predict_theta, h_value, shift, subgroup_size,
                          max_rl, block) {
  AEWMA <- 0
  D_lag1 <- 0
  D_lag2 <- 0

  for (i in seq_len(max_rl)) {
    x <- stats::rnorm(subgroup_size, mean = 0, sd = sqrt(shift))
    Mx <- stats::var(x)
    D <- .AEWMA_A + .AEWMA_B * log(Mx + .AEWMA_C)

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

    theta <- predict_theta(D, xt_lag1, xt_lag2)
    AEWMA <- theta * D + (1 - theta) * AEWMA
    D_lag2 <- D_lag1
    D_lag1 <- D

    if (abs(AEWMA) > h_value) return(i)
  }
  max_rl
}

.run_chains <- function(n_chains, lagged, predict_theta, h_value, shift,
                        subgroup_size, max_rl, block) {
  one <- if (lagged) .chain_lagged else .chain_blocked
  rl <- numeric(n_chains)
  for (j in seq_len(n_chains)) {
    rl[j] <- one(predict_theta, h_value, shift, subgroup_size, max_rl, block)
  }
  rl
}

# Executed once per stream, on a worker or in this process.
.stream_task <- function(i, streams, n_i, lagged, predict_theta, h_value, shift,
                         subgroup_size, max_rl, block) {
  assign(".Random.seed", streams[[i]], envir = globalenv())
  .run_chains(n_i[i], lagged, predict_theta, h_value, shift, subgroup_size,
              max_rl, block)
}

.make_streams <- function(seed, k) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(seed)
  streams <- vector("list", k)
  streams[[1L]] <- get(".Random.seed", envir = globalenv())
  for (i in seq_len(k - 1L)) {
    streams[[i + 1L]] <- parallel::nextRNGStream(streams[[i]])
  }
  streams
}

.simulate_run_lengths <- function(predict_theta, h_value, shift, subgroup_size,
                                  nsim, max_rl, seed, n_streams, n_cores,
                                  block = 500, lagged = FALSE) {
  .check_count(nsim, "nsim", min = 2)
  .check_count(max_rl, "max_rl")
  .check_count(block, "block")
  .check_count(n_streams, "n_streams")
  .check_count(n_cores, "n_cores")
  force(predict_theta)

  rng_state <- .save_rng()
  on.exit(.restore_rng(rng_state), add = TRUE)

  if (n_streams <= 1L) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(seed)
    return(.run_chains(nsim, lagged, predict_theta, h_value, shift,
                       subgroup_size, max_rl, block))
  }

  k <- as.integer(n_streams)
  streams <- .make_streams(seed, k)
  n_i <- rep(nsim %/% k, k) + (seq_len(k) <= nsim %% k)
  workers <- min(.effective_cores(n_cores), k)

  if (workers <= 1L) {
    rl <- lapply(
      seq_len(k), .stream_task, streams = streams, n_i = n_i, lagged = lagged,
      predict_theta = predict_theta, h_value = h_value, shift = shift,
      subgroup_size = subgroup_size, max_rl = max_rl, block = block
    )
    return(unlist(rl, use.names = FALSE))
  }

  old_threads <- .set_single_thread()
  on.exit(.restore_threads(old_threads), add = TRUE)
  cl <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  # workers must see the same library as this session
  parallel::clusterCall(cl, function(lib, pkgs) {
    .libPaths(lib)
    for (p in pkgs) loadNamespace(p)
    invisible(NULL)
  }, .libPaths(), c("mlAEWMA", "nnet", "e1071", "rpart", "gbm", "kernlab", "qrnn"))

  rl <- parallel::parLapply(
    cl, seq_len(k), .stream_task, streams = streams, n_i = n_i,
    lagged = lagged, predict_theta = predict_theta, h_value = h_value,
    shift = shift, subgroup_size = subgroup_size, max_rl = max_rl,
    block = block
  )
  unlist(rl, use.names = FALSE)
}
