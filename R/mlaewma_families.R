#' AEWMA control charts by model family
#'
#' Run the machine learning based adaptive EWMA chart for process dispersion
#' on subgroup data. Each subgroup variance is transformed as in Castagliola
#' (2005), a pre-fitted model predicts the smoothing weight `theta`
#' (restricted to \[0, 1\]) and the chart statistic is updated as
#' `AEWMA_t = theta_t * D_t + (1 - theta_t) * AEWMA_(t-1)`. The first
#' subgroup with `|AEWMA_t| > h` is the first signal.
#'
#' The four functions share all arguments. `mlaewma_neural()` offers `"nn"`,
#' `"qrnn"` and `"mcqrnn"`; `mlaewma_svm()` offers `"svr"`, `"svqr"` and
#' `"rsvqr"`; `mlaewma_tree()` offers `"dt"` and `"gb"`; `mlaewma_dqr()` has the
#' single model `"dqr"`, whose smoothing weight also depends on the two
#' previous values of `D`. With `method = "all"` every model of the family is
#' run and compared (for the one-model DQR family `"all"` is the same as
#' `"dqr"`). The `method` is not case sensitive. The models, the
#' transformation constants and the control limits are calibrated for
#' subgroups of size 5 only.
#'
#' `first_signal` is the number of the first subgroup that signals, counted
#' from the start of `data`; it is `NA` if the chart never signals.
#'
#' **Choosing `sigma0`.** `sigma0` must be the standard deviation of the
#' process while it is in control, estimated from a clean reference sample.
#' If it is left `NULL`, the standard deviation of all supplied observations is
#' used, which is inflated when `data` already contains an out-of-control
#' stretch and can hide the shift.
#'
#' **Limitations of the models.** The charts are calibrated for an in-control
#' ARL of about 370 (within roughly 2 percent; the bisection tolerance is
#' +/- 5). Their ability to detect changes differs a lot between models and
#' is not guaranteed:
#' * `"qrnn"` and `"mcqrnn"` assign almost no weight to large positive `D`;
#'   their average run length for a moderate variance increase (variance
#'   multiplier 1.5) is longer than the in-control one. They are not suitable
#'   for monitoring variance increases.
#' * `"nn"`, `"dt"`, `"gb"` and `"dqr"` detect variance increases but
#'   practically not decreases. `"svr"`, `"svqr"` and `"rsvqr"` react to both
#'   directions.
#' * The models were trained on values of the transformed statistic between 0
#'   and 1.9, but about half of the in-control values of `D` are negative, so
#'   the predicted weights are partly extrapolated.
#'
#' See the vignette for simulated comparisons.
#'
#' Simulation (needed when `find_h = TRUE` or `shifts` is given) is controlled
#' by `n_streams` and `n_cores`; see [calculate_arl()]. The results depend on
#' `seed`, `n_streams` and the number of simulations, not on `n_cores`.
#'
#' @param data Numeric vector of observations (cut into subgroups of
#'   `subgroup_size`; an incomplete last subgroup is dropped) or a matrix or
#'   data frame with one subgroup per row.
#' @param method Model to use, or `"all"` for every model of the family.
#' @param h Control limit. `NULL` uses the calibrated default. For
#'   `method = "all"` a vector with one value per model, named
#'   (case-insensitive) or in model order.
#' @param sigma0 In-control standard deviation of the process. If `NULL`, the
#'   standard deviation of all observations is used.
#' @param subgroup_size Subgroup size; only 5 is supported.
#' @param find_h If `TRUE`, search `h` by bisection so that the in-control ARL
#'   equals `target_arl0` (slow).
#' @param target_arl0 Target in-control ARL for `find_h`.
#' @param h_search_nsim Number of simulations per bisection step.
#' @param h_validation_nsim Number of simulations for the final validation of
#'   the found `h`.
#' @param arl_tolerance Bisection stops when the simulated ARL0 is within this
#'   distance of `target_arl0`.
#' @param shifts Optional vector of variance multipliers for which the
#'   out-of-control ARL is simulated.
#' @param shift_nsim,shift_max_rl Number of simulations and maximum run length
#'   used for `shifts`.
#' @param nsim,max_rl Simulation settings used by `find_h` and `shifts`
#'   (`nsim` is not used for the application run).
#' @param seed Random seed for simulations.
#' @param n_cores Number of worker processes used for simulation (default 1;
#'   at most 2 while `R CMD check` runs).
#' @param n_streams Number of independent random number streams used for
#'   simulation; the default 3 reproduces the calibration of the default
#'   control limits.
#' @return For a single model a list with `method`, `h`, `first_signal`,
#'   `results` (data frame with `Subgroup`, `Variance`, `D`, `Theta`, `AEWMA`,
#'   `h`, `Exceeds_Limit`), `calibration` (only with `find_h = TRUE`) and, if
#'   `shifts` is given, `performance`. For `method = "all"` a list with one such
#'   element per model plus `comparison`, `performance`, `best`,
#'   `rank_summary` and `calibration`.
#' @references Castagliola, P. (2005). A new S2-EWMA control chart for
#'   monitoring the process variance. *Quality and Reliability Engineering
#'   International*, 21(8), 781-794. \doi{10.1002/qre.686}
#'
#'   Noor-ul-Amin, M., Kazmi, M. W., Alkhalaf, S., Abdel-Khalek, S. and
#'   Nabi, M. (2024). Machine learning based parameter-free adaptive EWMA
#'   control chart to monitor process dispersion. *Scientific Reports*, 14,
#'   31271. \doi{10.1038/s41598-024-82699-8}
#'
#'   Kazmi, M. W. and Noor-ul-Amin, M. (2024). Adaptive EWMA control chart by
#'   using support vector regression. *Quality and Reliability Engineering
#'   International*, 40(7), 3831-3843. \doi{10.1002/qre.3603}
#'
#'   Capizzi, G. and Masarotto, G. (2003). An adaptive exponentially weighted
#'   moving average control chart. *Technometrics*, 45(3), 199-207.
#'   \doi{10.1198/004017003000000023}
#'
#'   Atar, G., Purutçuoğlu, V. and Bursalı, A. (2026). Comparison of machine
#'   learning models in a parameter-free adaptive EWMA control chart for
#'   process dispersion monitoring. Accepted at the 3rd International
#'   Congress on Data Science and Information Technologies (INFTEC 2026).
#' @examples
#' data(patient)
#' res <- mlaewma_neural(patient, method = "nn", sigma0 = 4.05)
#' res$first_signal
#' head(res$results)
#' mlaewma_tree(patient, method = "all", sigma0 = 4.05)$comparison
#'
#' # a variance increase after subgroup 20 (in-control sigma0 = 2)
#' set.seed(10)
#' x <- c(rnorm(100, sd = 2), rnorm(60, sd = 4))
#' mlaewma_svm(x, method = "svqr", sigma0 = 2)$first_signal
#' @name mlaewma_family
NULL

#' @rdname mlaewma_family
#' @export
mlaewma_neural <- function(data, method = "qrnn", h = NULL, sigma0 = NULL,
                         subgroup_size = 5, find_h = FALSE,
                         target_arl0 = 370, h_search_nsim = 3000,
                         h_validation_nsim = 20000, arl_tolerance = 5,
                         shifts = NULL, shift_nsim = 10000,
                         shift_max_rl = 10000, nsim = 500, max_rl = 10000,
                         seed = 42, n_cores = 1L, n_streams = 3L) {
  .aewma_family("neural", data, method, h, sigma0, subgroup_size, find_h,
                target_arl0, h_search_nsim, h_validation_nsim, arl_tolerance,
                shifts, shift_nsim, shift_max_rl, nsim, max_rl, seed, n_cores,
                n_streams)
}

#' @rdname mlaewma_family
#' @export
mlaewma_svm <- function(data, method = "svqr", h = NULL, sigma0 = NULL,
                      subgroup_size = 5, find_h = FALSE,
                      target_arl0 = 370, h_search_nsim = 3000,
                      h_validation_nsim = 20000, arl_tolerance = 5,
                      shifts = NULL, shift_nsim = 10000,
                      shift_max_rl = 10000, nsim = 500, max_rl = 10000,
                      seed = 42, n_cores = 1L, n_streams = 3L) {
  .aewma_family("svm", data, method, h, sigma0, subgroup_size, find_h,
                target_arl0, h_search_nsim, h_validation_nsim, arl_tolerance,
                shifts, shift_nsim, shift_max_rl, nsim, max_rl, seed, n_cores,
                n_streams)
}

#' @rdname mlaewma_family
#' @export
mlaewma_tree <- function(data, method = "gb", h = NULL, sigma0 = NULL,
                       subgroup_size = 5, find_h = FALSE,
                       target_arl0 = 370, h_search_nsim = 3000,
                       h_validation_nsim = 20000, arl_tolerance = 5,
                       shifts = NULL, shift_nsim = 10000,
                       shift_max_rl = 10000, nsim = 500, max_rl = 10000,
                       seed = 42, n_cores = 1L, n_streams = 3L) {
  .aewma_family("tree", data, method, h, sigma0, subgroup_size, find_h,
                target_arl0, h_search_nsim, h_validation_nsim, arl_tolerance,
                shifts, shift_nsim, shift_max_rl, nsim, max_rl, seed, n_cores,
                n_streams)
}

#' @rdname mlaewma_family
#' @export
mlaewma_dqr <- function(data, method = "dqr", h = NULL, sigma0 = NULL,
                      subgroup_size = 5, find_h = FALSE,
                      target_arl0 = 370, h_search_nsim = 3000,
                      h_validation_nsim = 20000, arl_tolerance = 5,
                      shifts = NULL, shift_nsim = 10000,
                      shift_max_rl = 10000, nsim = 500, max_rl = 10000,
                      seed = 42, n_cores = 1L, n_streams = 3L) {
  .aewma_family("dqr", data, method, h, sigma0, subgroup_size, find_h,
                target_arl0, h_search_nsim, h_validation_nsim, arl_tolerance,
                shifts, shift_nsim, shift_max_rl, nsim, max_rl, seed, n_cores,
                n_streams)
}
