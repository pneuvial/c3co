library("c3co")

context("Toy data generation")


test_that("output of getToyData has correct dimensions and contents", {
    len <- 100L   ## Number of loci
    K <- 2L       ## Number of subclones
    J <- 6L       ## Number of segments
    n <- 4L       ## Number of samples
    eps <- 1.0

    configs <- expand.grid(
        intercept = c(TRUE, FALSE),
        returnLocus = c(TRUE, FALSE),
        dimension = 1:2)
    for (cc in 1:nrow(configs)) {
        config <- configs[cc, ]
        intercept <- config[["intercept"]]
        returnLocus <- config[["returnLocus"]]
        dimension <- config[["dimension"]]
        dat <- getToyData(n, len, nbClones=K, nbSegs=J, eps = 0.0, 
                          dimension = dimension,
                          intercept = intercept,
                          returnLocus = returnLocus) 
        
        ## weights
        expect_equal(nrow(dat$W), n)
        expect_equal(ncol(dat$W), K)
        
        ## segment-level data
        sdat <- dat$segment
        if (dimension == 1) {
            expect_equal(nrow(sdat$Y), n)
            expect_equal(ncol(sdat$Y), J)
            expect_equal(nrow(sdat$Z), K)
            expect_equal(ncol(sdat$Z), J)
        } else {
            expect_equal(nrow(sdat$Y[[1]]), n)
            expect_equal(ncol(sdat$Y[[1]]), J)
            expect_equal(nrow(sdat$Z[[1]]), K)
            expect_equal(ncol(sdat$Z[[1]]), J)
            
            expect_equal(nrow(sdat$Y[[2]]), n)
            expect_equal(ncol(sdat$Y[[2]]), J)
            expect_equal(nrow(sdat$Z[[2]]), K)
            expect_equal(ncol(sdat$Z[[2]]), J)
        }    
        ## locus-level data
        ldat <- dat$locus
        if (returnLocus) {
            if (dimension == 1) {
                expect_equal(nrow(ldat$Y), n)
                expect_equal(ncol(ldat$Y), len)
                expect_equal(nrow(ldat$Z), K)
                expect_equal(ncol(ldat$Z), len)
            } else {
                expect_equal(nrow(ldat$Y[[1]]), n)
                expect_equal(ncol(ldat$Y[[1]]), len)
                expect_equal(nrow(ldat$Z[[1]]), K)
                expect_equal(ncol(ldat$Z[[1]]), len)

                expect_equal(nrow(ldat$Y[[2]]), n)
                expect_equal(ncol(ldat$Y[[2]]), len)
                expect_equal(nrow(ldat$Z[[2]]), K)
                expect_equal(ncol(ldat$Z[[2]]), len)
            }
        } else {
            expect_null(ldat)
        }        
    }    
})


test_that("getToyData throws errors as appropriate", {
    # less samples than clones:
    expect_error(getToyData(n = 5L, len = 10L, nbClones = 6L, nbSegs = 5L, eps = 0.0))
    # less segments than clones:
    expect_error(getToyData(n = 10L, len = 10L, nbClones = 5L, nbSegs = 3L, eps = 0.0))
})