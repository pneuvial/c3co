#' Function to load Facets data and transform them to c3co format
#'
#' @export
#' @param pathFacets The path to load Facets data csv.gz format.
#' @return A data frame under PSCBS format
loadFacetsdata <- function(pathFacets){
    ### Load Facets data
    dat <- lapply(list.files(pathFacets), FUN=function(ff) {
        df <- facets::readSnpMatrix(file.path(pathFacets, ff))
        xx=facets::preProcSample(df)
        dat=xx$pmat
        ## Rename chromosome, x and CT to segment with c3co
        df <- data.frame(chr=dat$chrom, pos=dat$maploc)
        df$tcn <- 2*dat$rCountT/dat$rCountN
        df$dh <- 2*abs(dat$vafT-1/2)
        df$dh[dat$het==0] <- NA
        return(df)
    })
    return(dat)
}

#' Function to transform Facets data, perfom the segmentation and save it into output.dir
#'
#' @export
#' @param pathFacets The path to load Facets data.
#' @param output.dir Directory to save segmentation
#' @param stat "TCN or "C1C2" paramater to segment the data. If \code{stat==TCN}, the segmentation will be done on TCN only. 
#' @return A list which contains the breakpoints by chromosome and also the binning of TCN, C1 and C2.
#' @examples
#' if (require("facets", quietly=TRUE)) {
#' pathFacets <- system.file("extdata",package="facets")
#' output.dir <- R.utils::Arguments$getWritablePath("output")
#' Facetswrapper(pathFacets,output.dir=output.dir, stat="TCN")
#' \dontrun{
#' resc3co <- c3co(NULL, pathSeg=output.dir)
#' }
#' }
Facetswrapper <- function (pathFacets, output.dir, stat){
    if (!requireNamespace("facets", quietly=TRUE)) {
        stop("Package 'facets' needed. Please install it from github/mskcc",             call. = FALSE)
    }
    ### To do may be cut this function into several function
    dat <- loadFacetsdata(pathFacets)
    ### Joint segmentation of all samples
    resSeg <-  segmentData(dat, stat=stat)
    ### Perform the c3co method
    saveRDS(resSeg, file=file.path(output.dir, "segDat.rds"))
    message(sprintf("segment data has been saved to %s in segDat.rds file\n",output.dir)) 
}
