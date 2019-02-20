context("Consistency of the optimization algorithm: get.W")

reference_data <- readRDS(system.file("extdata", "data4ConsistencyTests.rds", package = "c3co"))

test_that("Consistency of get.W", {

  n <- nrow(reference_data$W)
  K <- ncol(reference_data$W)
  J <- ncol(reference_data$Y)
  
  ## solving with c3co (limSolve)
  W_c3co <- c3co:::get.W(t(reference_data$Z), reference_data$Y)  
    
  expect_lt(sum((reference_data$W_hat_ref - W_c3co)^2), 1e-10)
})

