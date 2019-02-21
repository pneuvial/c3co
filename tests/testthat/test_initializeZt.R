context("Initialization of the latent features")

set.seed(7)
dataAnnotTP <- acnr::loadCnRegionData(dataSet = "GSE11976", tumorFraction = 1.0)
dataAnnotN  <- acnr::loadCnRegionData(dataSet = "GSE11976", tumorFraction = 0.0)

len <- 500 * 10  ## Number of loci
K <- 3L          ## Number of subclones
n <- 15L         ## Number of samples

bkps <- list(
  c(100, 250) * 10,
  c(150, 400) * 10,
  c(150, 400) * 10
)

regions <- list(
  c("(0,3)", "(0,2)", "(1,2)"),
  c("(1,1)", "(0,1)", "(1,1)"),
  c("(0,2)", "(0,1)", "(1,1)")
)

datSubClone <- buildSubclones(len = len,
                              nbClones = K,
                              bkps = bkps,
                              regions = regions,
                              dataAnnotTP = dataAnnotTP,
                              dataAnnotN = dataAnnotN)
stopifnot(is.list(datSubClone), length(datSubClone) == K)

W <- rSparseWeightMatrix(nb.samp = n, nb.arch = K)
stopifnot(identical(dim(W), c(n, K)))

dat <- mixSubclones(subClones = datSubClone, W = W)
stopifnot(is.list(dat), length(dat) == nrow(W))

seg <- segmentData(dat)
Y1 <- t(seg$Y1)
Y2 <- t(seg$Y2)
J  <- ncol(Y1)

test_that("Outputs of initializeZt have the expected dimensions", {
    
    ## flavors <- c("hclust", "nmf", "archetypes", "svd", "subsampling")
    flavors <- c("hclust", "nmf", "svd", "subsampling")
    
    for (ff in flavors) {
        test_that("output of initializeZt() has correct dimensions and contents", {
            Zt <- initializeZt(Y1, Y2, K = K, flavor = ff)
            expect_equal(nrow(Zt$Z1), J)
            expect_equal(nrow(Zt$Z2), J)
            
            expect_equal(ncol(Zt$Z1), K)
            expect_equal(ncol(Zt$Z2), K)
            
            Zt <- initializeZt(Y1, K=K, flavor = ff)
            expect_null(Zt$Z2)
            ##        expect_equal(Zt$Z1, Zt$Z)  ## There is no 'Z' /HB 2018-02-27
            
            expect_equal(nrow(Zt$Z1), J)
            expect_equal(ncol(Zt$Z1), K)
        })
    }
})
