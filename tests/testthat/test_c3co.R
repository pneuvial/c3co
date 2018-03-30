context("c3co test")

set.seed(7)
dataAnnotTP <- acnr::loadCnRegionData(dataSet = "GSE11976", tumorFraction = 1.0)
dataAnnotN <- acnr::loadCnRegionData(dataSet = "GSE11976", tumorFraction = 0.0)

len <- 500 * 10
nbClones <- 3

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
                              nbClones = nbClones,
                              bkps = bkps,
                              regions = regions,
                              dataAnnotTP = dataAnnotTP,
                              dataAnnotN = dataAnnotN)
stopifnot(is.list(datSubClone), length(datSubClone) == 3L)

M <- rSparseWeightMatrix(nb.samp = 15L, nb.arch = 3L)
stopifnot(all.equal(dim(M), c(15, 3)))

dat <- mixSubclones(subClones = datSubClone, W = M)
stopifnot(is.list(dat), length(dat) == nrow(M))

l1 <- 10^(-seq(from = 2, to = 8, by = 1))
parameters.grid <- list(lambda = l1, nb.arch = 2:6)

test_that("c3co terminates on C1,C2", {
  res <- c3co(dat, parameters.grid = parameters.grid)
  stopifnot(inherits(res, "c3coFit"),
            is.list(res@segDat),
            all(c("Y1", "Y2", "Y") %in% names(res@segDat)))
              
  df <- createZdf(res, chromosomes = 1L, idxBest = 3L)
  stopifnot(is.data.frame(df))
            
  pvePlot(res)
  Zplot(df)
})

test_that("c3co terminates on TCN", {
  resC <- c3co(dat, parameters.grid = parameters.grid, stat = "TCN")
  stopifnot(inherits(resC, "c3coFit"))
            
  df <- createZdf(resC, chromosomes = 1L, idxBest = 3L)
  stopifnot(is.data.frame(df))
  
  pvePlot(resC)
  Zplot(df)
})
