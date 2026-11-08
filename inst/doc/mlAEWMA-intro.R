## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(mlAEWMA)

## ----single-------------------------------------------------------------------
data(patient)
res <- mlaewma_neural(patient, method = "nn", sigma0 = 4.05)
res$h
res$first_signal
head(res$results, 3)

## ----signal-------------------------------------------------------------------
set.seed(10)
x <- c(rnorm(100, sd = 2), rnorm(60, sd = 4))
mlaewma_svm(x, method = "svqr", sigma0 = 2)$first_signal
mlaewma_dqr(x, sigma0 = 2)$first_signal

## ----compare------------------------------------------------------------------
cmp <- mlaewma_compare(patient, sigma0 = 4.05)
cmp$comparison

## ----arl----------------------------------------------------------------------
calculate_arl("nn", h_value = 0.846, nsim = 200, seed = 1)

## ----shifts-------------------------------------------------------------------
h <- c(nn = 0.846, qrnn = 1.968, svqr = 0.568)
sapply(c(1, 1.5, 2), function(s) {
  sapply(names(h), function(m) {
    round(calculate_arl(m, h[[m]], shift = s, nsim = 300, seed = 3)[["ARL"]])
  })
})

