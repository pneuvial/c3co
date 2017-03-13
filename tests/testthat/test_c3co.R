dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
len <- 500*10
nbClones <- 3
bkps <- list(c(100,250)*10, c(150,400)*10,c(150,400)*10)
regions <-list(c("(0,3)", "(0,2)","(1,2)"), c("(1,1)", "(0,1)","(1,1)"), c("(0,2)", "(0,1)","(1,1)"))
datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN, nbClones, bkps, regions)
M <- getWeightMatrix(100,0, 3, 15, sparse.coeff=0.7, contam.coeff=0.6, contam.max=2)
dat <- mixSubclones(subClones=datSubClone, M)
l1 <- seq(from=1e-6, to=1e-5, length=3)
l2 <- seq(from=1e-6, to=1e-5, length=3)

test_that("c3co terminates on C1,C2", {
    casResC1C2 <- c3co(dat, lambda1.grid=l1, lambda2.grid=l2, nb.arch.grid=2:6)
})

test_that("c3co terminates on TCN", {
    casRes <- c3co(dat, stat="TCN", lambda1.grid=l1, lambda2.grid=l2, nb.arch.grid=2:6)
})
