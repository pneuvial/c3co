library("c3co")

context("Initialization of the latent features")

dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
len <- 500*10   ## Number of loci
K <- 3L         ## Number of subclones
n <- 12L        ## Number of samples
bkps <- list(c(100, 250)*10, c(150, 400)*10, c(150, 400)*10)
stopifnot(length(bkps) == K)
J <- length(unique(unlist(bkps))) + 1L
regions <- list(c("(0,3)", "(0,1)", "(1,2)"),
                c("(1,1)", "(0,1)", "(1,1)"),
                c("(0,2)", "(0,1)", "(1,1)"))
datSubClone <- buildSubclones(len, nbClones=K, bkps=bkps, regions=regions,
                              dataAnnotTP=dataAnnotTP, dataAnnotN=dataAnnotN)
W <- rSparseWeightMatrix(nb.samp=n, nb.arch=K, sparse.coeff=0.90)
simu <- mixSubclones(subClones=datSubClone, W=W)
seg <- segmentData(simu)
Y1 <- t(seg$Y1)
Y2 <- t(seg$Y2)
expect_equal(dim(Y1), dim(Y2))
expect_equal(nrow(Y1), n)
expect_equal(ncol(Y1), J)

## flavors <- c("hclust", "nmf", "archetypes", "svd", "subsampling")
flavors <- c("hclust", "nmf", "svd", "subsampling")

for (ff in flavors) {
    test_that("output of initializeZt() has correct dimensions and contents", {
        Zt <- initializeZt(Y1, Y2, K=K, flavor=ff)
        expect_equal(nrow(Zt$Z1), J)
        expect_equal(nrow(Zt$Z2), J)
        
        expect_equal(ncol(Zt$Z1), K)
        expect_equal(ncol(Zt$Z2), K)
        
        Zt <- initializeZt(Y1, K=K, flavor=ff)
        expect_null(Zt$Z2)
##        expect_equal(Zt$Z1, Zt$Z)  ## There is no 'Z' /HB 2018-02-27
        
        expect_equal(nrow(Zt$Z1), J)
        expect_equal(ncol(Zt$Z1), K)
    })
}
