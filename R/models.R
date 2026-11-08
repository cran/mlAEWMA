# Prediction of the smoothing weight theta from the transformed statistic D.
# The fitted models live in R/sysdata.rda (object `.aewma_models`), built by
# data-raw/build_sysdata.R.

#' @importFrom nnet nnet
#' @importFrom e1071 svm
#' @importFrom rpart rpart
#' @importFrom gbm gbm
#' @importFrom kernlab ksvm kernelMult rbfdot
#' @importFrom qrnn qrnn.predict mcqrnn.predict
#' @importFrom stats predict var sd rnorm setNames aggregate
NULL

# Quantile SVM with a Gaussian kernel (formerly qrsvm:::predict.qrsvm).
# f(x) = sum_i alpha_i * exp(-sigma * ||x - x_i||^2) + b0
.predict_rsvqr <- function(model, D) {
  kernlab::kernelMult(
    kernlab::rbfdot(sigma = model$sigma),
    matrix(D, ncol = 1L),
    model$xtrain,
    model$alpha
  ) + model$b0
}

# theta must stay inside [0, 1].
.clamp_theta <- function(theta) {
  pmin(pmax(abs(as.numeric(theta)), 0), 1)
}

# Vectorised theta predictor for the eight memoryless models.
.theta_fun <- function(method) {
  model <- .aewma_models[[method]]
  raw <- switch(
    method,
    nn = function(D) stats::predict(model, data.frame(xt = D)),
    qrnn = function(D) qrnn::qrnn.predict(matrix(D, ncol = 1L), parms = model),
    mcqrnn = function(D) {
      # trained at tau = 0.1, 0.5, 0.9; the median is used as theta
      qrnn::mcqrnn.predict(matrix(D, ncol = 1L), parms = model, tau = 0.5)
    },
    svr = function(D) stats::predict(model, data.frame(xt = D)),
    svqr = function(D) kernlab::predict(model, data.frame(xt = D)),
    rsvqr = function(D) .predict_rsvqr(model, D),
    dt = function(D) stats::predict(model, data.frame(xt = D)),
    gb = function(D) stats::predict(model, data.frame(xt = D), n.trees = 500),
    stop("unknown method: ", method, call. = FALSE)
  )
  function(D) .clamp_theta(raw(D))
}

# DQR: theta depends on D and the two previous D values.
.theta_fun_dqr <- function() {
  cf <- .aewma_models[["dqr"]]
  function(D, xt_lag1, xt_lag2) {
    .clamp_theta(
      cf[["(Intercept)"]] + cf[["xt"]] * D +
        cf[["xt_lag1"]] * xt_lag1 + cf[["xt_lag2"]] * xt_lag2
    )
  }
}
