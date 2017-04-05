#' Joint segmentation
#'
#' Joint segmentation by Recursive Binary Segmentation followed by
#' Dynamic Programming
#'
#' @param dat A list of \code{data.frame} containing
#'  \describe{
#'   \item{tcn}{Total copy number}
#'   \item{dh}{Mirrored B allele fraction}
#'   \item{pos}{Position on the genome}
#'   \item{chr}{Chromosome}
#' }
#' 
#' @param stat "TCN or "C1C2" parameter to segment the data.
#' If \code{stat == TCN}, the segmentation will be done on TCN only.
#' 
#' @param verbose A logical value indicating whether to print extra
#' information. Defaults to FALSE
#' 
#' @return Binned Minor and Major copy number with list of breakpoints
#'
#' 
#' @references Gey, S., & Lebarbier, E. (2008). Using CART to Detect Multiple
#'   Change Points in the Mean for Large Sample.
#'   http://hal.archives-ouvertes.fr/hal-00327146/

#' @examples
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 3
#' bkps <- list(c(100, 250)*10, c(150, 400)*10, c(150, 400)*10)
#' regions <- list(c("(0,3)", "(0,2)", "(1,2)"), 
#' c("(1,1)", "(0,1)", "(1,1)"), c("(0,2)", "(0,1)", "(1,1)"))
#' datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN,
#'                               nbClones, bkps, regions)
#' M <- getWeightMatrix(100, 0, 3, 15, sparse.coeff=0.7,
#'                      contam.coeff=0.6, contam.max=2)
#' dat <- mixSubclones(subClones=datSubClone, M)
#' res <- segmentData(dat)
#' res2 <- segmentData(dat, stat="TCN")
#'
#' @importFrom matrixStats binMeans
#' @importFrom jointseg jointSeg
#' @export
segmentData <- function(dat, stat=c("C1C2", "TCN"), verbose=FALSE) {
    stat <- match.arg(stat)
    
    checkColNames <- lapply(dat, FUN=function(dd) {
        coln <- colnames(dd)
        ecn <- c("tcn", "dh", "pos", "chr") ## expected
        mm <- match(ecn, coln)
        if (any(is.na(mm))) {
            str <- sprintf("('%s')", paste(ecn, collapse="', '"))
            stop("Argument 'data' should contain columns named ", str)
        }
    })
    
    tcn <- lapply(dat, FUN = function(x) x$tcn)
    tcn <- Reduce(cbind, tcn)
    if (stat == "C1C2") {
        dh <- lapply(dat, FUN = function(x) x$dh)
        dh <- Reduce(cbind, dh)
        dataToSeg <- cbind(tcn, dh)
    } else if (stat == "TCN") {
        dataToSeg <- cbind(tcn) 
    } 
    chrs <- unique(dat[[1]]$chr)
    bkpPosByCHR <- list()
    Y1 <- Y2 <- NULL
    Y <- DH <- NULL
    for (cc in chrs) {
        if (verbose) {
            message(sprintf("chr %s", cc))
        }
        ww <- which(dat[[1]]$chr == cc)
        if (verbose) {
            message("Joint segmentation")
        }
        resSeg <- jointSeg(Y=dataToSeg[ww, ], method="RBS", K=100,
                           modelSelectionMethod="Birge")
        bkp <- resSeg$bestBkp
        pos <- dat[[1]]$pos[ww]
        bkpPos <- rowMeans(cbind(pos[bkp], pos[bkp+1]))
        xOut <- c(min(pos), bkpPos, max(pos))
        xOut <- sort(unique(xOut))
        
        tcn <- as.matrix(tcn)
        
        if (stat == "C1C2") {
            dh <- as.matrix(dh)
            binDatTCN <- matrix(NA_real_, nrow=length(xOut)-1, ncol=ncol(tcn))
            for (bb in seq_len(ncol(tcn))) {
                means <- binMeans(y=tcn[ww, bb], x=pos, bx=xOut, na.rm=TRUE)
                binDatTCN[, bb] <- means
            }
            binDatDH <- matrix(NA_real_, nrow=length(xOut)-1, ncol=ncol(dh))
            for (bb in seq_len(ncol(dh))) {
                means <- binMeans(y=dh[ww, bb], x=pos, bx=xOut, na.rm=TRUE)
                binDatDH[, bb] <- means
            }
            
            idxNAC1 <- which(rowSums(is.na(binDatTCN)) > 0)
            idxNAC2 <- which(rowSums(is.na(binDatDH)) > 0)
            idxNA <- unique(c(idxNAC1, idxNAC2))
            if (length(idxNA)) {
                binDatTCNwithoutNA <- binDatTCN[-idxNA, ]
                binDatDHwithoutNA <- binDatDH[-idxNA, ]
                bkpPosByCHR[[cc]] <- bkpPos[-idxNA]
            } else {
                binDatTCNwithoutNA <- binDatTCN
                binDatDHwithoutNA <- binDatDH
                bkpPosByCHR[[cc]] <- bkpPos
            }
            ## Define Y1 and Y2
            Y <- rbind(Y, binDatTCNwithoutNA)
            DH <- rbind(DH, binDatDHwithoutNA)
            Y1 <- Y*(1-DH)/2
            Y2 <- Y*(1+DH)/2
        } else {
            binDatTCN <- matrix(NA_real_, nrow=length(xOut)-1, ncol=ncol(tcn))
            for (bb in seq_len(ncol(tcn))) {
                means <- binMeans(y=tcn[ww, bb], x=pos, bx=xOut, na.rm=TRUE)
                binDatTCN[, bb] <- means
            }
            idxNA <- which(rowSums(is.na(binDatTCN)) > 0)
            if (length(idxNA)) {
                binDatTCNwithoutNA <- binDatTCN[-idxNA, ]
                bkpPosByCHR[[cc]] <- bkpPos[-idxNA]
            } else {
                bkpPosByCHR[[cc]] <- bkpPos
                binDatTCNwithoutNA <- binDatTCN
            }
            ## Define Y1 and Y2
            Y <- rbind(Y, binDatTCNwithoutNA)
            Y1 <- NA_real_
            Y2 <- NA_real_
        }
    }
    return(list(Y1=Y1, Y2=Y2, Y=Y, bkp=bkpPosByCHR))
}
