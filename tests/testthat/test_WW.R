set.seed(90)

dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)

len <- 500*10
nbClones <- 3
bkps <- list(c(100, 250)*10, c(150, 400)*10, c(150, 400)*10)
regions <- list(c("(0,3)", "(0,2)", "(1,2)"),
                c("(1,1)", "(0,1)", "(1,1)"), 
                c("(0,2)", "(0,1)", "(1,1)"))

datSubClone <- buildSubclones(len, nbClones, bkps, regions, dataAnnotTP, dataAnnotN)
M <- rSparseWeightMatrix(15, 3, sparse.coeff=0.7)
dat <- mixSubclones(subClones=datSubClone, M)
seg <- segmentData(dat)

test_that("fitC3co terminates on C1C2", {
  lambda <- 0.01
  Y <- list(Y1=t(seg$Y1), Y2=t(seg$Y2))
  Z0 <- initializeZ(Y$Y1, Y$Y2, p=3)
  Z <- Z0[c("Z1", "Z2")]
  res <- positiveFusedLasso(Y, Z, lambda=rep(lambda, 2), verbose=TRUE)
})

