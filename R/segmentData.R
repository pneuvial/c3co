#' Joint segmentation
#'
#' Joint segmentation by Recursive Binary Segmentation followed by Dynamic Programming
#'
#' @export
#' @param dat A list of \code{data.frame} containing
#'  \describe{
#'   \item{tcn}{Total copy number}
#'   \item{dh}{Mirrored B allele fraction}
#'   \item{pos}{Position on the genome}
#'   \item{chr}{Chromosome}
#'   }
#' @param stat "TCN or "C1C2" paramater to segment the data. If \code{stat==TCN}, the segmentation will be done on TCN only. 
#' @param verbose A logical value indicating whether to print extra information. Defaults to FALSE
#' @return Binned Minor and Major copy number with list of breakpoints
#' #' 
#' @references Gey, S., & Lebarbier, E. (2008). Using CART to Detect Multiple
#'   Change Points in the Mean for Large Sample.
#'   http://hal.archives-ouvertes.fr/hal-00327146/

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
#' dat <- mixSubclones(subClones=datSubClone, M)
#' res <- segmentData(dat)
#' res2 <- segmentData(dat, stat="TCN")
segmentData <- function(dat, stat=c("C1C2", "TCN"), verbose=FALSE, ...){
    stat <- match.arg(stat)
    
    checkColNames <- lapply(dat, function(dd) {
        coln <- colnames(dd)
        ecn <- c("tcn", "dh", "pos", "chr") ## expected
        mm <- match(ecn, coln)
        if (any(is.na(mm))) {
            str <- sprintf("('%s')", paste(ecn, collapse="','"))
            stop("Argument 'data' should contain columns named ", str)
        }
    })

    ## purrr::map_df seems to require a *named* list
    if (is.null(names(dat))) {
        names(dat) <- seq_along(dat)
    }
    tcn <- purrr::map_df(dat, function(x) x$tcn)
    if (stat=="C1C2") {
        dh <- purrr::map_df(dat, function(x) x$dh)
        dataToSeg <- cbind(tcn, dh)
    } else if (stat=="TCN") {
        dataToSeg <- tcn
    } 

    len <- nrow(dataToSeg)
    chrs <- unique(dat[[1]]$chr)
    bkpPosByCHR <- list()
    Y1 <- Y2 <- NULL
    Y <- NULL
    for (cc in chrs) {
        if (verbose) {
            message(sprintf("chr %s", cc))
        }
        ww <- which(dat[[1]]$chr==cc)
        if (verbose) {
            message("Joint segmentation")
        }
        resSeg <- jointseg::jointSeg(Y=dataToSeg[ww,], method="RBS", K=100, modelSelectionMethod="Birge")
        bkp <- resSeg$bestBkp
        pos <- dat[[1]]$pos[ww]
        bkpPos <-rowMeans(cbind(pos[bkp], pos[bkp+1]))
        start <- c(min(pos), sort(c(pos[bkp+1])))
        end <- c(sort(pos[bkp]), max(pos))
        xOut <- c(min(pos), bkpPos, max(pos))
        xOut <- sort(unique(xOut))

        tcn <- as.matrix(tcn)
        dh <- as.matrix(dh)

        if (stat=="C1C2") {
            datC1 <- tcn*(1-dh)/2
            datC2 <- tcn*(1+dh)/2

            binDatC1 <- matrix(NA_real_, nrow=length(xOut)-1, ncol=ncol(datC1))
            for (bb in 1:ncol(datC1)) {
                means <- matrixStats::binMeans(y=datC1[ww, bb], x=pos, bx=xOut, na.rm=TRUE)
                binDatC1[, bb] <- means
            }
            binDatC2 <- matrix(NA_real_, nrow=length(xOut)-1, ncol=ncol(datC2))
            for (bb in 1:ncol(datC2)) {
                means <- matrixStats::binMeans(y=datC2[ww, bb], x=pos, bx=xOut, na.rm=TRUE)
                binDatC2[, bb] <- means
            }

            idxNAC1 <- which(rowSums(is.na(binDatC1))>0)
            idxNAC2 <- which(rowSums(is.na(binDatC2))>0)
            idxNA <- unique(c(idxNAC1, idxNAC2))
            if (length(idxNA)) {
                binDatC1withoutNA <- binDatC1[-idxNA,]
                binDatC2withoutNA <- binDatC2[-idxNA,]
                bkpPosByCHR[[cc]] <- bkpPos[-idxNA]
            } else {
                binDatC1withoutNA <- binDatC1
                binDatC2withoutNA <- binDatC2
                bkpPosByCHR[[cc]] <- bkpPos
            }
            ## Define Y1 and Y2
            Y1 <- rbind(Y1, binDatC1withoutNA)
            Y2 <- rbind(Y2, binDatC2withoutNA)
            Y <- Y1+Y2
        } else {
            binDatC <- matrix(NA_real_, nrow=length(xOut)-1, ncol=ncol(tcn))
            for (bb in 1:ncol(tcn)) {
                means <- matrixStats::binMeans(y=tcn[ww, bb], x=pos, bx=xOut, na.rm=TRUE)
                binDatC[, bb] <- means
            }
            idxNA <- which(rowSums(is.na(binDatC))>0)
            if (length(idxNA)) {
                binDatCwithoutNA <- binDatC[-idxNA,]
                bkpPosByCHR[[cc]] <- bkpPos[-idxNA]
            } else {
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
