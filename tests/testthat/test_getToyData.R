library("c3co")

context("Toy data generation")


test_that("output of getToyData has correct dimensions and contents", {
    len <- 100
    nbClones <- 2
    nbBkps <- 5
    eps <- 1
    n <- 4
    dat <- getToyData(n, len, nbClones, nbBkps, eps=0)  ## noiseless
    
    ## locus-level data
    ldat <- dat$locus
    expect_equal(nrow(ldat$Y), n)
    expect_equal(ncol(ldat$Y), len)
    expect_equal(nrow(ldat$W), n)
    expect_equal(ncol(ldat$W), nbClones)
    expect_equal(nrow(ldat$Z), nbClones)
    expect_equal(ncol(ldat$Z), len)
    
    ## segment-level data
    sdat <- dat$segment
    expect_equal(nrow(sdat$Y), n)
    expect_equal(ncol(sdat$Y), nbBkps+1)
    expect_equal(nrow(sdat$W), n)
    expect_equal(ncol(sdat$W), nbClones)
    expect_equal(nrow(sdat$Z), nbClones)
    expect_equal(ncol(sdat$Z), nbBkps+1)
})