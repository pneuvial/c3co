context("Consistency of the initialization step")

reference_data <- readRDS(system.file("extdata", "data4ConsistencyTests.rds", package = "c3co"))

test_that("Consistency of get.W", {

  n <- nrow(reference_data$W)
  K <- ncol(reference_data$W)
  J <- ncol(reference_data$Y)
  Y <- reference_data$Y
  Z <- reference_data$Z
  W <- reference_data$W
  Ybar <- colMeans(Y)

  Z0_hat_uncentered <- initializeZt(Y, K = 4)$Z1
  Z0_hat_centered   <- initializeZt(scale(Y, Ybar, FALSE) , K = 4)$Z1

  error_model      <- sum((Y - tcrossprod(W,t(Z)))^2)
  error_uncentered <- sum((Y - tcrossprod(W,Z0_hat_uncentered) )^2)
  error_centered   <- sum((Y - scale(tcrossprod(W,Z0_hat_centered), -Ybar, FALSE) )^2)

  expect_equal(error_uncentered, error_centered)
  expect_gt(error_uncentered   , error_model)
  expect_gt(error_centered     , error_model)
})

