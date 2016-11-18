#' Function to load PSCBS data and transform them to InCaSCN format
#'
#' @export
#' @param pathSegPSCBS The path to load PSCBS data.
#' @return A data frame under PSCBS format
loadPSCBSdata <- function(pathSegPSCBS){
### Load PSCBS data
  dat <- lapply(list.files(pathSegPSCBS), function (ff) {
    df <- readRDS(file.path(pathSegPSCBS, ff))$data
    ## Rename chromosome, x and CT to segment with InCaSCN
    df$chr <- df$chromosome
    df$pos <- df$x
    df$tcn <- df$CT
    df$dh <- df$rho
    return(df)
  })
### Check positions
  chr <- unique(dat[[1]]$chr)
  posFull <- lapply(chr, function (cc) {
    pp <- Reduce(intersect, lapply(dat, function(dd){
      d <- subset(dd, chr==cc)
      d$pos
    }))
  })
### Reduce data
  dat <- lapply(dat, function (ff) {
    df <- do.call(rbind, lapply(chr, function (cc){
      d <- subset(ff, chr==cc)
      d <- subset(d, pos%in%posFull[[cc]])
    }))
    return(df)
  })
  return(dat)
}

#' Function to segment data
#'
#' @export
#' @param dat The path to load PSCBS data.
#' @param stat "TCN or "C1C2" paramater to segment the data. If \code{stat==TCN}, the segmentation will be done on TCN only. 
#' @return A list which contains the breakpoints by chromosome and also the binning of TCN, C1 and C2.
segmentThroughInCaSCN <- function(dat, stat){
### 2 options joint segmentation of TCN profiles or TCN+dh$
  if(stat=="TCN"){
    resSeg <- InCaSCN:::segmentData(dat, stat=stat)
### Compute the DoH by segments
    seg.rho <- do.call(cbind, lapply(dat, function (df){
      chr.grid <- unique(df$chr)
      Y <- do.call(c, lapply(chr.grid, function (cc){
        dfCHR <- subset(df, chr==cc)
        xOut <- c(min(dfCHR$pos), resSeg$bkp[[cc]], max(dfCHR$pos))
        means <- matrixStats::binMeans(y = dfCHR$rho, x = dfCHR$pos, bx = xOut, na.rm = TRUE)
      }))
      return(Y)
    }))
### Compute the minor and major copy by segments
    seg.Y1 <- resSeg$Y*(1-seg.rho)/2
    seg.Y2 <- resSeg$Y*(1+seg.rho)/2
    
    resSeg$Y1 <- seg.Y1
    resSeg$Y2 <- seg.Y2
  }else if(stat=="C1C2"){
    ## This solution seems to be better (more breakpoints)
    resSeg <- InCaSCN:::segmentData(dat, stat=stat)
  }else{
    stop("stat must be TCN or C1C2")
  }
  return(resSeg)
}

#' Function to transform PSCBS data, perfom the segmentation and save it into output.dir
#'
#' @export
#' @param pathSegPSCBS The path to load PSCBS data.
#' @param output.dir Directory to save segmentation
#' @param stat "TCN or "C1C2" paramater to segment the data. If \code{stat==TCN}, the segmentation will be done on TCN only. 
#' @return A list which contains the breakpoints by chromosome and also the binning of TCN, C1 and C2.
PSCBSwrapper <- function (pathSegPSCBS,output.dir, stat){
  dat <- loadPSCBSdata(pathSegPSCBS)
### Joint segmentation of all samples
  resSeg <- segmentThroughInCaSCN(dat, stat)
  saveRDS(resSeg, file.path(output.dir, "segDat.rds"))
  cat(sprintf("segment data has been saved to %s in segDat.rds file\n",output.dir)) 
}
