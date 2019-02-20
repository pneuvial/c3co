context("Consistency of the optimization algorithm")

reference_data <- readRDS(system.file("extdata", "data4ConsistencyTests.rds", package = "c3co"))

test_that("Consistency of get.Z", {
  
  n <- nrow(reference_data$W)
  K <- ncol(reference_data$W)
  J <- ncol(reference_data$Y)
  
  ## solving with c3co (glmnet)
  WtWm1 <- tcrossprod(backsolve(qr.R(qr(reference_data$W)), x = diag(K)))
  ## FIXME: include the penalty factor in get.Zt
  Z_c3co <- t(c3co:::get.Zt(Y, lambda=lambda/(2*n*J), W = reference_data$W, WtWm1 = WtWm1))

  expect_lt(sum((reference_data$Zt_hat_ref - Z_c3co)^2), 1e-4)
})

