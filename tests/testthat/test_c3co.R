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

M <- getWeightMatrix(prop.max = 100, prob.min = 0,
                     nb.arch = 3L, nb.samp = 15L,
                     sparse.coeff = 0.7,
                     contam.coeff = 0.6,
                     contam.max = 2)

dat <- mixSubclones(subClones = datSubClone, W = M)

l1 <- seq(from = 1e-6, to = 1e-5, length.out = 3L)
l2 <- seq(from = 1e-6, to = 1e-5, length.out = 3L)

parameters.grid <- list(lambda1 = l1, lambda2 = l2, nb.arch = 2:6)

test_that("c3co terminates on C1,C2", {
  res <- c3co(dat, parameters.grid = parameters.grid)
  df <- createZdf(res, chromosomes=1, idxBest=2)
  Zplot(df)
})

test_that("c3co terminates on TCN", {
  resC <- c3co(dat, parameters.grid = parameters.grid, stat = "TCN")
  df <- createZdf(resC, chromosomes=1, idxBest=2)
  Zplot(df)
})
