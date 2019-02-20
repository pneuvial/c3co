library("c3co")

context("get.W")

test_that("get.W recovers W in noiseless situations", {
    configs <- expand.grid(
        sigSize = 100,
        sigDim = 1:2,
        nbClones = 1:3,
        nbSegs = c(1, 2, 8),
        nbSamples = c(1, 2, 6)
    )
    
    for (cc in 1:nrow(configs)) {
        config <- configs[cc, ]
        n <- config[["nbSamples"]]
        K <- config[["nbClones"]]
        J <- config[["nbSegs"]]
        M <- config[["sigDim"]]
        if ((n < K) || (J < K)) {
            expect_error(getToyData(n = n, len = config[["sigSize"]], 
                                    nbClones = K, nbSegs = J, 
                                    dimension = M, eps = 0))
        } else {
            dat <- getToyData(n = n, len = config[["sigSize"]], 
                              nbClones = K, nbSegs = J, 
                              dimension = M, eps = 0)
            W <- dat$W
            Y <- dat$loc$Y
            Z <- dat$loc$Z
            
            if (M > 1) {
                ## stack dimensions (if any)
                Y <- do.call(cbind, args = Y)
                Z <- do.call(cbind, args = Z)
            }    
            Zt <- t(Z)
            What <- c3co:::get.W(Zt, Y)
            expect_lt(max((W - What)^2), 1e-2)
        }
    }
})

