# Reference values were produced by the original, unpackaged scripts
# before the package was created.
D <- ref("ref_D.rds")
patient <- D$patient
s <- D$sigma0

test_that("application results equal the reference for every family", {
  ap <- ref("ref_app_patient.rds")
  expect_equal(mlaewma_neural(patient, "all", sigma0 = s), ap$neural, tolerance = 1e-12)
  expect_equal(mlaewma_svm(patient, "all", sigma0 = s), ap$svm, tolerance = 1e-12)
  expect_equal(mlaewma_tree(patient, "all", sigma0 = s), ap$tree, tolerance = 1e-12)
  expect_equal(mlaewma_dqr(patient, "dqr", sigma0 = s), ap$dqr, tolerance = 1e-12)
})

test_that("D is scale invariant", {
  aps <- ref("ref_app_patient_scaled.rds")
  expect_equal(mlaewma_neural(patient / s, "all", sigma0 = 1), aps$neural, tolerance = 1e-12)
  r1 <- mlaewma_neural(patient, "nn", sigma0 = s)$results$D
  r2 <- mlaewma_neural(patient / s, "nn", sigma0 = 1)$results$D
  expect_lt(max(abs(r1 - r2)), 1e-10)
})

test_that("theta stays in [0, 1] and equals the reference grid", {
  tg <- ref("ref_theta_grid.rds")
  for (m in setdiff(names(tg$theta), "dqr")) {
    th <- mlAEWMA:::.theta_fun(m)(tg$grid)
    expect_equal(th, tg$theta[[m]], tolerance = 1e-12, info = m)
    expect_true(all(th >= 0 & th <= 1), info = m)
  }
  g <- tg$grid
  l1 <- c(0, g[-length(g)])
  l2 <- c(0, 0, g[seq_len(length(g) - 2)])
  expect_equal(mlAEWMA:::.theta_fun_dqr()(g, l1, l2), tg$theta$dqr, tolerance = 1e-12)
})

test_that("rsvqr predictor reproduces the former qrsvm predict method", {
  rs <- ref("ref_rsvqr_pred.rds")
  expect_lt(max(abs(mlAEWMA:::.predict_rsvqr(mlAEWMA:::.aewma_models$rsvqr, rs$grid) - rs$pred)), 1e-10)
})

test_that("mlaewma_compare equals the reference and keeps its structure", {
  cp <- mlaewma_compare(patient, sigma0 = s)
  expect_equal(cp, ref("ref_compare_patient.rds"), tolerance = 1e-12)
  expect_named(cp, c("comparison", "models", "performance", "rank_summary", "family_results"))
  expect_equal(nrow(cp$comparison), 9)
})

test_that("input is validated", {
  expect_error(mlaewma_neural(patient, "nn", sigma0 = -1), "sigma0")
  expect_error(mlaewma_neural(patient[, 1:4], "nn", sigma0 = s), "4 column")
  expect_error(mlaewma_neural(c(1, NA, 3, 4, 5), "nn", sigma0 = 1), "missing")
  expect_error(mlaewma_neural(1:3, "nn", sigma0 = 1), "complete subgroup")
  expect_error(mlaewma_neural(patient, "svr", sigma0 = s))
  expect_error(mlaewma_compare(patient, h = c(nn = 1)), "9 names")
  expect_error(calculate_arl("nn", 0.8, nsim = 10, n_cores = 0), "n_cores")
  expect_error(calculate_arl("nn", 0.8, nsim = 1.5), "nsim")
})

test_that("odd inputs give clear errors instead of NaN or hangs", {
  expect_error(mlaewma_neural(rep(5, 50)), "no variation")
  expect_error(mlaewma_neural(patient * 1e200, "nn", sigma0 = 2e200), "not finite")
  expect_error(mlaewma_neural(patient, "nn", sigma0 = s, h = NA), "`h`")
  expect_error(mlaewma_neural(patient, "nn", sigma0 = s, h = -1), "`h`")
  expect_error(mlaewma_neural(patient, "nn", sigma0 = s, h = c(1, 2)), "`h`")
  expect_error(mlaewma_neural(patient, "all", sigma0 = s, h = c(1, NA, 1)), "`h`")
  expect_error(mlaewma_neural(patient, "nn", sigma0 = s, shifts = c(1, NA)), "shifts")
  expect_error(mlaewma_neural(patient, "nn", sigma0 = s, shifts = -1), "shifts")
  expect_error(mlaewma_neural(patient, "nn", sigma0 = s, seed = NA), "seed")
  expect_error(mlaewma_neural(patient, "nn", sigma0 = s, find_h = NA), "find_h")
  expect_error(mlaewma_neural(factor(rep(1:5, 4)), "nn", sigma0 = 1), "numeric")
  expect_error(mlaewma_neural(matrix(1:100, 5, 20), "nn", sigma0 = 1), "20 column")
  expect_error(mlaewma_neural(patient, "gb", sigma0 = s), "mlaewma_neural")
  expect_error(calculate_arl("nn", 0.8, nsim = 10, block = 0), "block")
  expect_error(calculate_arl("nn", 0.8, nsim = 10, shift = NA), "shift")
  expect_error(calculate_arl("nn", NA, nsim = 10), "h_value")
  expect_error(calculate_arl("nn", 0.8, nsim = 10, seed = NULL), "seed")
  expect_error(calculate_arl("nn", 0.8, nsim = 1), "nsim")
})

test_that("method is case-insensitive, ts input and DQR 'all' work", {
  expect_equal(mlaewma_neural(patient, "NN", sigma0 = s)$h, 0.846)
  x <- as.numeric(t(patient))
  expect_equal(mlaewma_neural(stats::ts(x), "nn", sigma0 = s)$results,
               mlaewma_neural(patient, "nn", sigma0 = s)$results)
  expect_equal(mlaewma_dqr(patient, "all", sigma0 = s)$h, 2.742)
})

test_that("a variance increase is detected", {
  set.seed(10)
  x <- c(rnorm(100, sd = 2), rnorm(60, sd = 4))
  expect_equal(mlaewma_dqr(x, sigma0 = 2)$first_signal, 22)
  expect_true(mlaewma_svm(x, "svqr", sigma0 = 2)$first_signal > 20)
})

test_that("block does not change the results", {
  a <- calculate_arl("nn", 0.846, nsim = 40, seed = 3, block = 1)
  b <- calculate_arl("nn", 0.846, nsim = 40, seed = 3, block = 500)
  expect_identical(a, b)
})

test_that("user's RNG kind produces no warning", {
  old <- suppressWarnings(RNGkind(sample.kind = "Rounding"))
  on.exit(suppressWarnings(RNGkind(sample.kind = old[3])))
  expect_no_warning(calculate_arl("nn", 0.846, nsim = 10))
})
