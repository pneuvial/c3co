context("Consistency of the optimization algorithm: get.Zt")

reference_data <- readRDS(system.file("extdata", "data4ConsistencyTests.rds", package = "c3co"))

test_that("Consistency of get.Zt", {
  
  n <- nrow(reference_data$W)
  K <- ncol(reference_data$W)
  J <- ncol(reference_data$Y)
  
  ## solving with c3co (glmnet)
  WtWm1 <- tcrossprod(backsolve(qr.R(qr(reference_data$W)), x = diag(K)))
  Z_c3co <- t(c3co:::get.Zt(reference_data$Y, lambda = reference_data$lambda, W = reference_data$W, WtWm1 = WtWm1))

  expect_lt(sum((reference_data$Z_hat_ref - Z_c3co)^2), 1e-3)
})

