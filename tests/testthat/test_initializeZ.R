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
    test_that("output of initializeZ has correct dimensions and contents", {
        res <- initializeZ(Y1, Y2, p=p, flavor=ff)
        expect_equal(nS, nrow(res$Z1))
        expect_equal(nS, nrow(res$Z2))
        
        expect_equal(p, ncol(res$Z1))
        expect_equal(p, ncol(res$Z2))
        
        res <- initializeZ(Y1, p=nbClones, flavor=ff)
        expect_null(res$Z2)
        expect_equal(res$Z1, res$Z)
        
        expect_equal(nS, nrow(res$Z))
        expect_equal(p, ncol(res$Z))
    })
}