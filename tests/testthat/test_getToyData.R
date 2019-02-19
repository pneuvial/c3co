library("c3co")

context("Toy data generation")


test_that("output of getToyData has correct dimensions and contents", {
    len <- 100
    nbClones <- 2
    nbBkps <- 5
    eps <- 1
    n <- 4
    
    for (intercept in c(TRUE, FALSE)) {
        for (returnLocus in c(TRUE, FALSE)) {
            dat <- getToyData(n, len, nbClones, nbBkps, eps=0, 
                              intercept = intercept,
                              returnLocus = returnLocus) 
            
            ## weights
            expect_equal(nrow(dat$W), n)
            expect_equal(ncol(dat$W), nbClones + intercept)
            
            ## locus-level data
            ldat <- dat$locus
            if (returnLocus) {
                expect_equal(nrow(ldat$Y), n)
                expect_equal(ncol(ldat$Y), len)
                expect_equal(nrow(ldat$Z), nbClones + intercept)
                expect_equal(ncol(ldat$Z), len)
            } else {
                expect_null(ldat)
            }        
            ## segment-level data
            sdat <- dat$segment
            expect_equal(nrow(sdat$Y), n)
            expect_equal(ncol(sdat$Y), nbBkps+1)
            expect_equal(nrow(sdat$Z), nbClones + intercept)
            expect_equal(ncol(sdat$Z), nbBkps+1)
        }    
    }
})