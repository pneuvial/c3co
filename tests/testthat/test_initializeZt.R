library("c3co")

context("Initialization of the latent features")

dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
len <- 500*10
nbClones <- 3
bkps <- list(c(100, 250)*10, c(150, 400)*10, c(150, 400)*10)
regions <- list(c("(0,3)", "(0,1)", "(1,2)"),
                c("(1,1)", "(0,1)", "(1,1)"),
                c("(0,2)", "(0,1)", "(1,1)"))
datSubClone <- buildSubclones(len, nbClones, bkps, regions, dataAnnotTP, dataAnnotN)
M <- rSparseWeightMatrix(12, nbClones, 0.90)
simu <- mixSubclones(subClones=datSubClone, M)
seg <- segmentData(simu)
Y1 <- t(seg$Y1)
Y2 <- t(seg$Y2)

## flavors <- c("hclust", "nmf", "archetypes", "svd", "subsampling")
flavors <- c("hclust", "nmf", "svd", "subsampling")
nS <- ncol(Y1)
p <- 3

for (ff in flavors) {
    test_that("output of initializeZt() has correct dimensions and contents", {
        Zt <- initializeZt(Y1, Y2, p=p, flavor=ff)
        expect_equal(nS, nrow(Zt$Z1))
        expect_equal(nS, nrow(Zt$Z2))
        
        expect_equal(p, ncol(Zt$Z1))
        expect_equal(p, ncol(Zt$Z2))
        
        Zt <- initializeZt(Y1, p=nbClones, flavor=ff)
        expect_null(Zt$Z2)
##        expect_equal(Zt$Z1, Zt$Z)  ## There is no 'Z' /HB 2018-02-27
        
        expect_equal(nS, nrow(Zt$Z1))
        expect_equal(p, ncol(Zt$Z1))
    })
}
