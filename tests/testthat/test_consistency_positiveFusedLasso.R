context("Consistency of the optimization algorithm: postitiveFusedLasso")

reference_data <- readRDS(system.file("extdata", "data4ConsistencyTests.rds", package = "c3co"))

test_that("Consistency of positiveFusedLasso", {

  tol <- 1e-6

  ## if we stick to one iteration, we call get.Wt once so this should 
  ## match "reference_data$W_hat_ref" just like in test_consistency_getw
  pfl <- positiveFusedLasso(
    Y      = list(reference_data$Y)   ,
    Zt     = list(t(reference_data$Z)),
    lambda = reference_data$lambda, eps = 1e-1, max.iter = 1L, verbose = TRUE
  )
  What  <- slot(pfl, "W")
  expect_lt(sum((What - reference_data$W_hat_ref)^2), tol)

  ## Now we can call get.Zt with the current estimation What and check 
  ## that it is matching the value returned by postitiveFusedLasso
  K      <- ncol(reference_data$Z)
  WtWm1  <- tcrossprod(backsolve(qr.R(qr(What)), x = diag(K)))
  Zt_ref <- c3co:::get.Zt(reference_data$Y, reference_data$lambda, What, WtWm1)
  Zthat  <- slot(pfl, "Zt")$Z1

  expect_lt(sum((Zthat - Zt_ref)^2), tol)
})
