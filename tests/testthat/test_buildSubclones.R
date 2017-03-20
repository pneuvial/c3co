library("c3co")

context("Construction of subclones")

dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
len <- 500*10
nbClones <- 2
bkps <- list(c(100,250)*10, c(150,400)*10)
regions <-list(c("(0,3)", "(0,2)","(1,2)"), c("(1,1)", "(0,1)","(1,1)"))
datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN, nbClones, bkps, regions)
names(datSubClone) <- seq_along(datSubClone)

test_that("built subclones have the expected size and column names", {
    expect_length(datSubClone, nbClones)
    enames <- c("ct", "baft", "genotype", "region","cn","bafn", "pos")
    for (ii in seq_along(datSubClone)) {
        sc <- datSubClone[[ii]]
        expect_equal(nrow(sc), len)
        
        nms <- names(sc)
        expect_equal(ncol(sc), length(enames))
        for (enm in enames) {
            expect_true(enm %in% nms)
        }
    }
})


test_that("built subclones have identical genotypes", {
    geno <- datSubClone[[1]]$genotype
    for (ii in seq_along(datSubClone)) {
        gg <- datSubClone[[ii]]$genotype
        expect_identical(gg, geno)
    }
})