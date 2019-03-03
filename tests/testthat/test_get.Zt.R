library("c3co")

context("get.Z")

test_that("get.Z recovers Z in (almost) noiseless situations for small lambda", {
    ## randomness in getToyData sometimes gives results out of tolerance band
    # set.seed(2) ## Gives error in R-devel (R 3.6.0)
    set.seed(3)
    
    lambdas = c(0, 1e-5)
    eps <- 1e-8  ## avoids glmnet complaining about 0 variance at standardization
    tol <- 5e-3
    
    configs <- expand.grid(
        sigSize = 20,
        sigDim = 1:2,
        nbClones = 1:3,
        nbSegs = c(1, 2, 8),
        nbSamples = c(2, 6))
    
    for (cc in 1:nrow(configs)) {
        config <- configs[cc, ]
        n <- config[["nbSamples"]]
        K <- config[["nbClones"]]
        J <- config[["nbSegs"]]
        M <- config[["sigDim"]]
        
        if ((n < K) || (J < K)) {
            expect_error(getToyData(n = n, len = config[["sigSize"]], 
                                    nbClones = K, nbSegs = J, 
                                    dimension = M, eps = eps))
        } else {
            dat <- getToyData(n = n, len = config[["sigSize"]], 
                              nbClones = K, nbSegs = J, 
                              dimension = M, eps = eps)
            W <- dat$W
            Y <- dat$loc$Y
            Z <- dat$loc$Z

            WtWm1 <- tcrossprod(backsolve(qr.R(qr(W)), x=diag(K)))
            
            for (lambda in lambdas) {
                if (M == 1) {
                    Z1that <- c3co:::get.Zt(Y, lambda = lambda, W = W, WtWm1 = WtWm1)
                    Z1hat <- t(Z1that)
                    expect_lt(max((Z - Z1hat)^2), tol)
                } else {
                    expect_equal(M, 2)
                    Z1that <- c3co:::get.Zt(Y[[1]], lambda = lambda, W = W, WtWm1 = WtWm1)
                    Z1hat <- t(Z1that)
                    expect_lt(max((Z[[1]] - Z1hat)^2), tol)
                    
                    Z2that <- c3co:::get.Zt(Y[[2]], lambda = lambda, W = W, WtWm1 = WtWm1)
                    Z2hat <- t(Z2that)
                    expect_lt(max((Z[[2]] - Z2hat)^2), tol)
                }
            }
        }
    }
})


