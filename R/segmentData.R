#' Segmentation function
#'
#' @export
#' @param dat A list of data frame for each patient containing the total copy number \code{tcn}, the mirrored B allele fraction \code{dh} and \code{chr} and \code{pos}.
#' @param stat "TCN or "C1C2" paramater to segment the data. If \code{stat==TCN}, the segmentation will be done on TCN only. 
#' @return Binned Minor and Major copy number with list of breakpoints
#' @examples
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 3
#' bkps <- list(c(100,250)*10, c(150,400)*10,c(150,400)*10)
#' regions <-list(c("(0,3)", "(0,2)","(1,2)"), 
#' c("(1,1)", "(0,1)","(1,1)"), c("(0,2)", "(0,1)","(1,1)"))
#' datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN, nbClones, bkps, regions)
#' M <- getWeightMatrix(100,0, 3, 15, sparse.coeff=0.7, contam.coeff=0.6, contam.max=2)
#' dat <- apply(M, 1, mixSubclones, subClones=datSubClone, fracN=NULL)
#' res <- segmentData(dat)
#' res2 <- segmentData(dat, stat="TCN")
segmentData <- function(dat, stat="C1C2"){
  
  a <- lapply(dat, function (dd) {
    coln <- colnames(dd)
    ecn <- c("tcn", "dh", "pos", "chr") ## expected
    mm <- match(ecn, coln)
    if (any(is.na(mm))) {
      str <- sprintf("('%s')", paste(ecn, collapse="','"))
      stop("Argument 'data' should contain columns named ", str)
    }
  })
  
  YTCNtoSeg <- t(sapply(dat, function(cc) cc$tcn))
  YDHtoSeg <- t(sapply(dat, function(cc) cc$dh))

  if(stat=="C1C2"){
    dataToSeg <- t(rbind(YTCNtoSeg,YDHtoSeg))
  }else{
    dataToSeg <- t(YTCNtoSeg)
  }
  len <- nrow(dataToSeg)
  chrs <- unique(dat[[1]]$chr)
  bkpPosByCHR <- list()
  Y1 <- Y2 <- NULL
  Y <- NULL
  for(cc in chrs){
    message(sprintf("chr %s", cc))
    ww <- which(dat[[1]]$chr==cc)
### Segmentation step on TCN and DH
    message("segmentation step")
    resSeg <- jointseg::jointSeg(Y=dataToSeg[ww,],K=100, modelSelectionMethod="Birge")
    message("end segmentation")
    bkp <- resSeg$bestBkp
    pos <- dat[[1]]$pos[ww]
    bkpPos <-rowMeans(cbind(pos[bkp], pos[bkp+1]))
    start <- c(min(pos), sort(c(pos[bkp+1])))
    end <- c(sort(pos[bkp]), max(pos))
    xOut <- c(min(pos), bkpPos, max(pos))
    xOut <- sort(unique(xOut))

    if(stat=="C1C2"){
      datC1 <- YTCNtoSeg*(1-YDHtoSeg)/2
      datC2 <- YTCNtoSeg*(1+YDHtoSeg)/2

      binDatC1 <- matrix(NA_real_, nrow=length(xOut)-1, ncol=nrow(datC1))
      for (bb in 1:nrow(datC1)) {
        means <- matrixStats::binMeans(y=datC1[bb, ww], x=pos, bx=xOut, na.rm=TRUE)
        binDatC1[, bb] <- means
      }
      binDatC2 <- matrix(NA_real_, nrow=length(xOut)-1, ncol=nrow(datC2))
      for (bb in 1:nrow(datC2)) {
        means <- matrixStats::binMeans(y=datC2[bb, ww], x=pos, bx=xOut, na.rm=TRUE)
        binDatC2[, bb] <- means
      }

      idxNAC1 <- sort(unique(unlist(apply(cbind(binDatC1), 2, function(ll) which(is.na(ll))))))
      idxNAC2 <- sort(unique(unlist(apply(cbind(binDatC2), 2, function(ll) which(is.na(ll))))))
      idxNA <- unique(c(idxNAC1,idxNAC2))
      if(length(idxNA)){
        binDatC1withoutNA <- binDatC1[-idxNA,]
        binDatC2withoutNA <- binDatC2[-idxNA,]
        bkpPosByCHR[[cc]] <- bkpPos[-idxNA]
      }else{
        bkpPosByCHR[[cc]] <- bkpPos
        binDatC1withoutNA <- binDatC1
        binDatC2withoutNA <- binDatC2
      }
      ## Define Y1 and Y2
      Y1 <- rbind(Y1, binDatC1withoutNA)
      Y2 <- rbind(Y2, binDatC2withoutNA)
      Y <- Y1+Y2
    }else{
      binDatC <- matrix(NA_real_, nrow=length(xOut)-1, ncol=nrow(YTCNtoSeg))
      for (bb in 1:nrow(YTCNtoSeg)) {
        means <- matrixStats::binMeans(y=YTCNtoSeg[bb, ww], x=pos, bx=xOut, na.rm=TRUE)
        binDatC[, bb] <- means
      }
      idxNA <- sort(unique(unlist(apply(cbind(binDatC), 2, function(ll) which(is.na(ll))))))
      if(length(idxNA)){
        binDatCwithoutNA <- binDatC[-idxNA,]
        bkpPosByCHR[[cc]] <- bkpPos[-idxNA]
      }else{
        bkpPosByCHR[[cc]] <- bkpPos
        binDatCwithoutNA <- binDatC
      }
      ## Define Y1 and Y2
      Y <- rbind(Y, binDatC)
      Y1 <- NA
      Y2 <- NA
    }
  }
  return(list(Y1=Y1, Y2=Y2, Y=Y, bkp=bkpPosByCHR))
}
