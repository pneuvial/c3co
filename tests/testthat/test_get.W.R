library("c3co")

context("getZ")

test_that("get.W recovers W for a single normal clone", {
    dat <- getToyData(n = 10, len = 20, nbClones = 1, 
                      nbSegs = 2, dimension = 1, eps = 0)
    W <- dat$W
    Y <- dat$loc$Y
    Z <- dat$loc$Z
    Zt <- t(Z)
    What <- c3co:::get.W(Zt, Y)
    What
    
})

test_that("get.W recovers W in noiseless situations", {
    configs <- expand.grid(
        sigSize = 100,
        sigDim = 1:2,
        nbClones = 2:3,
        nbSegs = c(1, 2, 8),
        nbSamples = c(2, 4, 6)  ## cannot test with nbSamples = 1 if nbClones cannot be set to 1!
    )
    
    configs <- subset(configs, nbSamples >= nbClones)
    configs <- subset(configs, nbSegs >= nbClones)
    
    for (cc in 1:nrow(configs)) {
        config <- configs[cc, ]
        dime <- config[["sigDim"]]
        dat <- getToyData(n = config[["nbSamples"]], len = config[["sigSize"]], 
                          nbClones = config[["nbClones"]], nbSegs = config[["nbSegs"]], 
                          dimension = dime, eps = 0)
        W <- dat$W
        Y <- dat$loc$Y
        Z <- dat$loc$Z
        
        if (dime > 1) {
            ## stack dimensions (if any)
            Y <- do.call(cbind, args = Y)
            Z <- do.call(cbind, args = Z)
        }    
        Zt <- t(Z)
        What <- c3co:::get.W(Zt, Y)
        expect_lt(max((W - What)^2), 1e-2)
    }
})

