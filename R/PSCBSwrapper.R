#' Loads PSCBS data and transforms them into c3co format
#'
#' @param PSCBSdata A list that contains PSCBS data as produced by
#' the \pkg{PSCBS} package.
#'
#' @return A data frame under PSCBS format.
#'
#' @export
loadPSCBSdata <- function(PSCBSdata) {
    ### Load PSCBS data
    dat <- lapply(PSCBSdata, FUN=function(ff) {
        df <- ff$data
        ## Rename chromosome, x and CT to segment with c3co
        df$chr <- df$chromosome
        df$pos <- df$x
        df$tcn <- df$CT
        df$dh <- df$rho
        df
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
        as.data.frame(df)
    })
    dat
}

#' Transforms and aligns PSCBS data followed by a joint segmentation
#'
#' @param PSCBSdata A list that contains PSCBS data
#'
#' @param stat `"TCN"` or `"C1C2"` parameter to segment the data.
#' If `stat == "TCN"`, the segmentation will be done on TCN only.
#'
#' @param align If `TRUE`, the PSCBS list elements are expanded before
#' segmentation to have the exact same set of (chr, pos) loci, which is the
#' union of all elements. See [alignLoci()] for details.
#' 
#' @return A list which contains the breakpoints by chromosome and also the
#' binning of TCN, C1, and C2.
#'
#' @export
PSCBSwrapper <- function(PSCBSdata, stat, align = FALSE) {
  stopifnot(is.list(PSCBSdata))

  ## Joint segmentation of all samples
  dat <- loadPSCBSdata(PSCBSdata)

  ## Align?
  if (align) {
    dat <- alignLoci(dat)
  }

  ## Resegment
  resSeg <- segmentData(dat, stat = stat)
  
  ## Sanity checks
  stopifnot(ncol(resSeg$Y) == length(dat))
  
  resSeg
}
