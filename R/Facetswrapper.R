#' Function to load Facets data and transform them to InCaSCN format
#'
#' @export
#' @param pathFacets The path to load Facets data csv.gz format.
#' @return A data frame under PSCBS format
loadFacetsdata <- function(pathFacets){
### Load Facets data
  dat <- lapply(list.files(pathFacets), function (ff) {
    df <- facets::readSnpMatrix(file.path(pathFacets, ff))
    xx=facets::preProcSample(rcmat)
    dat=xx$pmat
    ## Rename chromosome, x and CT to segment with InCaSCN
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
Facetswrapper <- function (pathFacets,output.dir, stat){
### To do may be cut this function into several function
  dat <- loadFacetsdata(pathFacets)
### Joint segmentation of all samples
  resSeg <- segmentThroughInCaSCN(dat, output.dir, stat)
### Perform the InCaSCN method
  saveRDS(resSeg, file.path(output.dir, "segDat.rds"))
  cat(sprintf("segment data has been saved to %s in segDat.rds file\n",output.dir)) 
}
