#' Function to load PSCBS data and transform them to c3co format
#'
#' @param PSCBSdata A list that contains PSCBS data.
#'
#' @param pattern if you have several patients in the pathSegPSCBS directory.
#'
#' @return A data frame under PSCBS format
#'
#' @export
loadPSCBSdata <- function(PSCBSdata, pattern=NULL) {
    ### Load PSCBS data
    dat <- lapply(PSCBSdata, FUN=function(ff) {
        df <- ff$data
        ## Rename chromosome, x and CT to segment with c3co
        df$chr <- df$chromosome
        df$pos <- df$x
        df$tcn <- df$CT
        df$dh <- df$rho
        return(df)
    })
    ### Check positions
    chr <- unique(dat[[1]]$chr)
    posFull <- lapply(chr, FUN=function(cc) {
        pp <- Reduce(intersect, lapply(dat, FUN=function(dd) {
            d <- subset(dd, chr == cc)
            d$pos
        }))
    })
    ### Reduce data
    dat <- lapply(dat, FUN=function(ff) {
        df <- do.call(rbind, args=lapply(chr, function(cc) {
            d <- subset(ff, chr == cc)
            pos <- NULL; rm(list = "pos");
            d <- subset(d, pos %in% posFull[[cc]])
        }))
        return(as.data.frame(df))
    })
    return(dat)
}

#' Function to transform PSCBS data, perfom the segmentation and save it
#' into output.dir
#'
#' @param PSCBSdata A list that contains PSCBS data
#'
#' @param output.dir Directory to save segmentation
#'
#' @param stat "TCN or "C1C2" paramater to segment the data.
#' If \code{stat == TCN}, the segmentation will be done on TCN only.
#'
#' @param pattern if you have several patients in the pathSegPSCBS directory.
#'
#' @return A list which contains the breakpoints by chromosome and also the
#' binning of TCN, C1 and C2.
#'
#' @export
PSCBSwrapper <- function(PSCBSdata, pattern=NULL, output.dir, stat) {
    dat <- loadPSCBSdata(PSCBSdata, pattern)
    ### Joint segmentation of all samples
    resSeg <- segmentData(dat, stat=stat)
    saveRDS(resSeg, file=file.path(output.dir, "segDat.rds"))
    mprintf("segment data has been saved to %s in segDat.rds file", output.dir)
}
