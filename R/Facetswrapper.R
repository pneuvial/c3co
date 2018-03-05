#' Loads FACETS data and transforms them to c3co format
#'
#' @param pathFacets The path to a directory containing FACETS data files
#' (\file{*.csv.gz}) produced by the \pkg{facets} package.
#'
#' @return A data frame under PSCBS format.
#'
#' @export
loadFacetsdata <- function(pathFacets) {
    dat <- lapply(list.files(pathFacets), FUN=function(ff) {
        df <- facets::readSnpMatrix(file.path(pathFacets, ff))
        xx <- facets::preProcSample(df)
        dat <- xx$pmat
        ## Rename chromosome, x and CT to segment with c3co
        df <- data.frame(chr=dat$chrom, pos=dat$maploc)
        df$tcn <- 2*dat$rCountT/dat$rCountN
        df$dh <- 2*abs(dat$vafT-1/2)
        df$dh[dat$het == 0] <- NA_real_
        df
    })
    dat
}

#' Transforms FACETS data and performs a joint segmentation
#'
#' @param pathFacets The path to a directory containing FACETS data files
#' (\file{*.csv.gz}) produced by the \pkg{facets} package.
#'
#' @param stat `"TCN"` or `"C1C2"` paramater to segment the data.
#' If `stat == "TCN"`, the segmentation will be done on TCN only.
#'
#' @return A list which contains the breakpoints by chromosome and also the
#' binning of TCN, C1, and C2.
#'
#' @examples
#' if (require("facets", quietly=TRUE)) {
#' pathFacets <- system.file("extdata", package="facets")
#' segDat <- Facetswrapper(pathFacets, stat="TCN")
#' \dontrun{
#' resc3co <- c3co(NULL, segDat=segDat)
#' }
#' }
#'
#' @export
Facetswrapper <- function(pathFacets, stat) {
    if (!requireNamespace("facets", quietly=TRUE)) {
      stop("Package 'facets' needed. See https://github.com/mskcc/facets",
           call. = FALSE)
    }
    ### To do may be cut this function into several function
    dat <- loadFacetsdata(pathFacets)
    ### Joint segmentation of all samples
    resSeg <- segmentData(dat, stat=stat)
    ## Sanity checks
    stopifnot(ncol(resSeg$Y) == length(dat))
    resSeg
}
