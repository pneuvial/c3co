dataAnnotTP <- acnr::loadCnRegionData(dataSet = "GSE11976", tumorFraction = 1.0)
dataAnnotN <- acnr::loadCnRegionData(dataSet = "GSE11976", tumorFraction = 0.0)

len <- 500 * 10
nbClones <- 3L

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
                              dataAnnotTP = dataAnnotTP,
                              dataAnnotN = dataAnnotN,
                              nbClones = nbClones,
                              bkps = bkps,
                              regions = regions)

M <- rSparseWeightMatrix(15, 3)

dat <- mixSubclones(subClones = datSubClone, W = M)

l1 <- 10^(-seq(from = 2, to = 8, by = 1))

parameters.grid <- list(lambda = l1, nb.arch = 2:6)

test_that("c3co terminates on C1,C2", {
  res <- c3co(dat, parameters.grid = parameters.grid)
  df <- createZdf(res, chromosomes=1, idxBest=3)
  pvePlot(res)
  Zplot(df)
})

test_that("c3co terminates on TCN", {
  resC <- c3co(dat, parameters.grid = parameters.grid, stat = "TCN")
  df <- createZdf(resC, chromosomes=1, idxBest=3)
  pvePlot(resC)
  Zplot(df)
})
