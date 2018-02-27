context("c3co test on PSCBS data")

if (requireNamespace("c3co.data")) {
  data("PSCBSdata", package = "c3co.data")

  test_that("Produce (C1,C2) segmentation from PSCBS locus-level data", {
    segDat <- PSCBSwrapper(PSCBSdata, stat = "C1C2")
    stopifnot(is.list(segDat),
              all(c("Y1", "Y2", "Y", "bkp") %in% names(segDat)))
  }
  
  set.seed(7)
  lambda.grid <- seq(from = 1e-4, to = 1e-3, length.out = 10)
  parameters.grid <- list(lambda = lambda.grid, nb.arch = 2:5)

  test_that("c3co terminates on PSCBS locus-level data", {
    res <- c3co(dat = NULL, parameters.grid, segDat = segDat)
  })
}
