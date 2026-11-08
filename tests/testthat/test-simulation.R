test_that("ARL0 simulation equals the original scripts (3 streams)", {
  ar <- ref("ref_arl0_small.rds")
  for (nm in grep("_nc3$", names(ar), value = TRUE)) {
    x <- ar[[nm]]
    method <- sub("_nc3$", "", nm)
    if (method == "rsvqr") skip_on_cran()
    got <- calculate_arl(method, x$h, nsim = x$nsim, seed = 42, n_streams = 3)
    expect_equal(got, x$value, tolerance = 1e-12, info = nm)
  }
})

test_that("single-stream simulation equals the original single-core runs", {
  ar <- ref("ref_arl0_small.rds")
  for (nm in c("nn_nc1", "dt_nc1", "dqr_nc1")) {
    x <- ar[[nm]]
    got <- calculate_arl(sub("_nc1$", "", nm), x$h, nsim = x$nsim, seed = 42,
                         n_streams = 1)
    expect_equal(got, x$value, tolerance = 1e-12, info = nm)
  }
})

test_that("results do not depend on n_cores and are repeatable", {
  a <- calculate_arl("nn", 0.846, nsim = 60, seed = 11)
  expect_identical(a, calculate_arl("nn", 0.846, nsim = 60, seed = 11))
  # n_cores > 1 needs the installed package on the workers
  skip_on_cran()
  skip_if_not(dir.exists(file.path(system.file(package = "mlAEWMA"), "Meta")),
              "package is not installed (devtools::load_all)")
  b <- calculate_arl("nn", 0.846, nsim = 60, seed = 11, n_cores = 2)
  expect_identical(a, b)
})

test_that("the session's random number generator is left untouched", {
  set.seed(7)
  before <- .Random.seed
  invisible(calculate_arl("nn", 0.846, nsim = 20))
  invisible(calculate_arl("dqr", 2.742, nsim = 20, n_streams = 1))
  expect_identical(.Random.seed, before)
})

test_that("shifts and find_h paths equal the reference", {
  skip_on_cran()
  D <- ref("ref_D.rds")
  patient <- D$patient
  s <- D$sigma0
  p <- ref("ref_paths.rds")
  sh <- c(1, 1.5)
  expect_equal(mlaewma_neural(patient, "all", sigma0 = s, shifts = sh, shift_nsim = 100),
               p$neural_shift, tolerance = 1e-12)
  expect_equal(mlaewma_dqr(patient, "dqr", sigma0 = s, shifts = sh, shift_nsim = 100),
               p$dqr_shift, tolerance = 1e-12)
  r <- suppressWarnings(suppressMessages(
    mlaewma_tree(patient, "gb", sigma0 = s, find_h = TRUE,
               h_search_nsim = 100, h_validation_nsim = 200)))
  expect_equal(r, p$gb_findh, tolerance = 1e-12)
})
