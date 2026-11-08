# Frozen constants. Do not change: the models, the calibrated control limits
# and the Castagliola (2005) transformation constants belong together.

.AEWMA_A <- -0.8969
.AEWMA_B <- 2.3647
.AEWMA_C <- 0.5979

# Control limits h calibrated to an in-control ARL of about 370 with
# seed = 42 and n_streams = 3 (see the data-raw scripts).
.aewma_default_h <- c(
  nn = 0.846, qrnn = 1.968, mcqrnn = 1.968,
  svr = 0.661, svqr = 0.568, rsvqr = 0.703,
  dt = 1.385, gb = 1.355,
  dqr = 2.742
)

.ALL_METHODS <- names(.aewma_default_h)

# Upper bracket of the bisection used when find_h = TRUE.
.aewma_h_upper <- c(
  nn = 6, qrnn = 6, mcqrnn = 15,
  svr = 6, svqr = 6, rsvqr = 6,
  dt = 6, gb = 6,
  dqr = 10
)

.aewma_families <- list(
  neural = c("nn", "qrnn", "mcqrnn"),
  svm = c("svr", "svqr", "rsvqr"),
  tree = c("dt", "gb"),
  dqr = "dqr"
)
