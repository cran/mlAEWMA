# Helpers that keep the user's session untouched.

# Remember the state of the random number generator of the calling session.
.save_rng <- function() {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  list(had = had, old = old, kind = RNGkind())
}

# Put the random number generator back exactly as .save_rng() found it.
.restore_rng <- function(state) {
  # a user's non-default sample.kind (e.g. "Rounding") would warn here
  suppressWarnings(RNGkind(state$kind[1L], state$kind[2L], state$kind[3L]))
  if (state$had) {
    assign(".Random.seed", state$old, envir = globalenv())
  } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  invisible(NULL)
}

# Limit BLAS/OpenMP threads for the duration of a simulation.
.thread_vars <- c(
  "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
)

.set_single_thread <- function() {
  old <- Sys.getenv(.thread_vars, unset = NA)
  do.call(Sys.setenv, as.list(stats::setNames(rep("1", length(.thread_vars)), .thread_vars)))
  old
}

.restore_threads <- function(old) {
  unset <- names(old)[is.na(old)]
  if (length(unset)) Sys.unsetenv(unset)
  keep <- old[!is.na(old)]
  if (length(keep)) do.call(Sys.setenv, as.list(keep))
  invisible(NULL)
}

# Number of worker processes actually started. Never more than two while
# R CMD check is running.
.effective_cores <- function(n_cores) {
  n <- as.integer(n_cores)
  lim <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (nzchar(lim) && lim != "false") n <- min(n, 2L)
  max(n, 1L)
}

.check_count <- function(x, name, min = 1) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < min ||
        x != round(x)) {
    stop("`", name, "` must be a single whole number >= ", min, ".",
         call. = FALSE)
  }
  invisible(as.integer(x))
}

# single finite number > 0
.check_pos <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("`", name, "` must be a single finite number > 0.", call. = FALSE)
  }
  invisible(x)
}

.check_seed <- function(seed) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed != round(seed) || abs(seed) > .Machine$integer.max) {
    stop("`seed` must be a single whole number.", call. = FALSE)
  }
  invisible(seed)
}

.check_shifts <- function(shifts) {
  if (!is.numeric(shifts) || length(shifts) < 1L || any(!is.finite(shifts)) ||
        any(shifts <= 0)) {
    stop("`shifts` must be a numeric vector of finite values > 0.",
         call. = FALSE)
  }
  invisible(shifts)
}
