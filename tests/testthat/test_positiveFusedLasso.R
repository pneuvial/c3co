library("c3co")

context("positiveFusedLasso")

test_that("positiveFusedLasso recovers the truth in (almost) noiseless situations for small lambda", {
    lambdas = c(0, 1e-5)
    eps <- 1e-8  ## avoids glmnet complaining about 0 variance at standardization
    tol <- 5e-3
    set.seed(3)  ## randomness in getToyData sometimes gives results out of tolerance band
    
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
            if (M == 1) {
                Y <- list(Y)
                Z <- list(Z)
            }
            Zt <- lapply(Z, t)
            
            for (lambda in lambdas) {
                pfl <- positiveFusedLasso(Y, Zt = Zt, lambda = rep(lambda, M), 
                                          eps = 1e-1, max.iter = 50L)
                What <- pfl@W
                Yhat <- pfl@E
                Zhat <- pfl@Zt
                
                expect_lt(max((What - W)^2), tol)
                for (mm in 1:M) {
                    expect_lt(max((Yhat[[mm]] - Y[[mm]])^2), tol)
                    expect_lt(max((Zhat[[mm]] - Zt[[mm]])^2), tol)
                }
            }
        }
    }
})


test_that("positiveFusedLasso with intercept recovers the truth in (almost) noiseless situations for small lambda", {
    lambdas = c(0, 1e-5)
    eps <- 1e-8  ## avoids glmnet complaining about 0 variance at standardization
    tol <- 5e-3
    set.seed(3)  ## randomness in getToyData sometimes gives results out of tolerance band
    
    configs <- expand.grid(
        sigSize = 20,
        sigDim = 1:2,
        nbClones = 2:3,
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
            if (M == 1) {
                Y <- list(Y)
                Z <- list(Z)
            }
            Zt <- lapply(Z, t)
            
            for (lambda in lambdas) {
                cat("\n", n, K, J, M, lambda)
                pfl <- positiveFusedLasso(Y, Zt = Zt, lambda = rep(lambda, M), 
                                          intercept = TRUE,
                                          eps = 1e-1, max.iter = 50L)
                What <- pfl@W
                Yhat <- pfl@E
                Zhat <- pfl@Zt

                expect_lt(max((What - W)^2), tol)
                for (mm in 1:M) {
                    expect_lt(max((Yhat[[mm]] - Y[[mm]])^2), tol)
                    expect_lt(max((Zhat[[mm]] - sweep(Zt[[mm]], 2, colMeans(Zt[[mm]])))^2), tol)
                }
            }
        }
    }
})

