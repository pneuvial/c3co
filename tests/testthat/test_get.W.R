library("c3co")

context("getZ")

test_that("get.W recovers W in noiseless situations", {
    configs <- expand.grid(
        sigSize = 100,
        sigDim = 1:2,
        nbClones = 1:3,
        nbSegs = c(1, 2, 8),
        nbSamples = c(2, 4, 6)  ## cannot test with nbSamples = 1 if nbClones cannot be set to 1!
    )
    
    for (cc in 1:nrow(configs)) {
        config <- configs[cc, ]
        dime <- config[["sigDim"]]
        if ((config[["nbSamples"]] < config[["nbClones"]]) || (config[["nbSegs"]] < config[["nbClones"]])) {
            expect_error(getToyData(n = config[["nbSamples"]], len = config[["sigSize"]], 
                                    nbClones = config[["nbClones"]], nbSegs = config[["nbSegs"]], 
                                    dimension = dime, eps = 0))
        } else {
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
    }
})

