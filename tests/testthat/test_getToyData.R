library("c3co")

context("Toy data generation")


test_that("output of getToyData has correct dimensions and contents", {
    len <- 100
    nbClones <- 2
    nbSegs <- 6
    eps <- 1
    n <- 4

    configs <- expand.grid(
        intercept = c(TRUE, FALSE),
        returnLocus = c(TRUE, FALSE),
        dimension = 1:2)
    for (cc in 1:nrow(configs)) {
        config <- configs[cc, ]
        intercept <- config[["intercept"]]
        returnLocus <- config[["returnLocus"]]
        dimension <- config[["dimension"]]
        dat <- getToyData(n, len, nbClones, nbSegs, eps = 0, 
                          dimension = dimension,
                          intercept = intercept,
                          returnLocus = returnLocus) 
        
        ## weights
        expect_equal(nrow(dat$W), n)
        expect_equal(ncol(dat$W), nbClones)
        
        ## segment-level data
        sdat <- dat$segment
        if (dimension == 1) {
            expect_equal(nrow(sdat$Y), n)
            expect_equal(ncol(sdat$Y), nbSegs)
            expect_equal(nrow(sdat$Z), nbClones)
            expect_equal(ncol(sdat$Z), nbSegs)
        } else {
            expect_equal(nrow(sdat$Y[[1]]), n)
            expect_equal(ncol(sdat$Y[[1]]), nbSegs)
            expect_equal(nrow(sdat$Z[[1]]), nbClones)
            expect_equal(ncol(sdat$Z[[1]]), nbSegs)
            
            expect_equal(nrow(sdat$Y[[2]]), n)
            expect_equal(ncol(sdat$Y[[2]]), nbSegs)
            expect_equal(nrow(sdat$Z[[2]]), nbClones)
            expect_equal(ncol(sdat$Z[[2]]), nbSegs)
        }    
        ## locus-level data
        ldat <- dat$locus
        if (returnLocus) {
            if (dimension == 1) {
                expect_equal(nrow(ldat$Y), n)
                expect_equal(ncol(ldat$Y), len)
                expect_equal(nrow(ldat$Z), nbClones)
                expect_equal(ncol(ldat$Z), len)
            } else {
                expect_equal(nrow(ldat$Y[[1]]), n)
                expect_equal(ncol(ldat$Y[[1]]), len)
                expect_equal(nrow(ldat$Z[[1]]), nbClones)
                expect_equal(ncol(ldat$Z[[1]]), len)

                expect_equal(nrow(ldat$Y[[2]]), n)
                expect_equal(ncol(ldat$Y[[2]]), len)
                expect_equal(nrow(ldat$Z[[2]]), nbClones)
                expect_equal(ncol(ldat$Z[[2]]), len)
            }
        } else {
            expect_null(ldat)
        }        
    }    
})


test_that("getToyData throws errors as appropriate", {
    # less samples than clones:
    expect_error(getToyData(n = 5, len = 10, nbClones = 6, nbSegs = 5, eps = 0))
    # less segments than clones:
    expect_error(getToyData(n = 10, len = 10, nbClones = 5, nbSegs = 3, eps = 0))
})