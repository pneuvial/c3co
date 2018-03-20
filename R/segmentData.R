#' Joint segmentation
#' 
#' Joint segmentation by Recursive Binary Segmentation followed by Dynamic 
#' Programming
#' 
#' @param dat A list of data frames containing columns:
#' \describe{
#'   \item{`tcn`}{Total copy number}
#'   \item{`dh`}{Mirrored B allele fraction}
#'   \item{`pos`}{Position on the genome}
#'   \item{`chr`}{Chromosome}
#' }
#'   
#' @param stat `"TCN"` or `"C1C2"` parameter to segment the data.
#' If `stat == `"TCN"`, the segmentation will be done on TCN only.
#'   
#' @param verbose A logical value indicating whether to print extra
#' information. Defaults to `FALSE`.
#'   
#' @return Binned minor and major copy numbers with list of breakpoints.
#'   
#' @references Gey, S., & Lebarbier, E. (2008). Using CART to Detect Multiple 
#'   Change Points in the Mean for Large Sample. 
#'   http://hal.archives-ouvertes.fr/hal-00327146/
#'   
#' @references Pierre-Jean, M., Rigaill, G. & Neuvial, P. (2015). Performance 
#'   evaluation of DNA copy number segmentation methods.  Briefings in 
#'   Bioinformatics 16 (4): 600-615
#'   
#' @details This function is a wrapper around the [jointseg::jointSeg()] of
#' the \pkg{jointseg} package.
#'   
#' @examples
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 3
#' bkps <- list(c(100, 250)*10, c(150, 400)*10, c(150, 400)*10)
#' regions <- list(c("(0,3)", "(0,2)", "(1,2)"),
#' c("(1,1)", "(0,1)", "(1,1)"), c("(0,2)", "(0,1)", "(1,1)"))
#' datSubClone <- buildSubclones(len, nbClones, bkps, regions, dataAnnotTP, dataAnnotN)
#' M <- rSparseWeightMatrix(15, 3, sparse.coeff=0.7)
#' dat <- mixSubclones(subClones=datSubClone, M)
#' res <- segmentData(dat)
#' res2 <- segmentData(dat, stat="TCN")
#'
#' @importFrom matrixStats binMeans
#' @importFrom matrixStats rowAnyNAs
#' @importFrom jointseg jointSeg
#' @export
segmentData <- function(dat, stat=c("C1C2", "TCN"), verbose=FALSE) {
    stat <- match.arg(stat)

    checkColNames <- lapply(dat, FUN=function(dd) {
        coln <- colnames(dd)
        ecn <- c("tcn", "dh", "pos", "chr") ## expected
        mm <- match(ecn, coln)
        if (anyNA(mm)) {
            stop("Argument 'data' should contain columns named ",
                 comma(sQuote(ecn)))
        }
    })

    chrs <- unique(dat[[1]]$chr)
    stopifnot(!anyNA(chrs))
    
    ## Assert that all samples are for the same set of loci, which is assumed below
    if (length(dat) > 1) {
      chr1 <- dat[[1]]$chr
      x1 <- dat[[1]]$x
      for (ii in 2:length(dat)) {
        chr <- dat[[ii]]$chr
        if (length(chr) != length(chr1)) {
          stop(sprintf("Sample #%d is for different number of loci than Sample #1: %d != %d",
                       ii, length(chr), length(chr1)))
        }
        if (!all(chr == chr1, na.rm = TRUE)) {
          stop(sprintf("Sample #%d is for a different set of chromosomes than Sample #1", ii))
        }
        x <- dat[[ii]]$x
        if (!all(x == x1, na.rm = TRUE)) {
          stop(sprintf("Sample #%d is for a different set of positions than Sample #1", ii))
        }
      }
    }

    tcn <- lapply(dat, FUN = function(x) x$tcn)
    tcn <- Reduce(cbind, tcn)
    tcn <- as.matrix(tcn)

    if (stat == "C1C2") {
        dh <- lapply(dat, FUN = function(x) x$dh)
        dh <- Reduce(cbind, dh)
        dh <- as.matrix(dh)
        dataToSeg <- cbind(tcn, dh)
    } else if (stat == "TCN") {
        dataToSeg <- cbind(tcn)
    }

    names <- names(dat)
    
    bkpPosByCHR <- list()
    Y1 <- Y2 <- NULL
    Y <- DH <- NULL
    for (cc in chrs) {
        if (verbose) {
            mprintf("chr %s\n", cc)
        }
        ww <- which(dat[[1]]$chr == cc)
        if (verbose) message("Joint segmentation")

        stopifnot(length(ww) > 0, length(ww) >= 3L)
        resSeg <- jointSeg(Y=dataToSeg[ww, ], method="RBS", K=100,
                           modelSelectionMethod="Birge")
        bkp <- resSeg$bestBkp
        pos <- dat[[1]]$pos[ww]
        bkpPos <- rowMeans(cbind(pos[bkp], pos[bkp+1]))  ## FIXME: work on segments to avoid arbitrary bkp position?
        xOut <- c(min(pos), bkpPos, max(pos))
        xOut <- sort(unique(xOut))

        if (stat == "C1C2") {
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

            idxNAC1 <- which(rowAnyNAs(binDatTCN))
            idxNAC2 <- which(rowAnyNAs(binDatDH))
            idxNA <- unique(c(idxNAC1, idxNAC2))
            if (length(idxNA) > 0) {
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
            colnames(Y) <- names
            DH <- rbind(DH, binDatDHwithoutNA)
            Y1 <- Y*(1-DH)/2
            Y2 <- Y*(1+DH)/2
        } else {
            binDatTCN <- matrix(NA_real_, nrow=length(xOut)-1, ncol=ncol(tcn))
            for (bb in seq_len(ncol(tcn))) {
                means <- binMeans(y=tcn[ww, bb], x=pos, bx=xOut, na.rm=TRUE)
                binDatTCN[, bb] <- means
            }
            idxNA <- which(rowAnyNAs(binDatTCN))
            if (length(idxNA) > 0) {
                binDatTCNwithoutNA <- binDatTCN[-idxNA, ]
                bkpPosByCHR[[cc]] <- bkpPos[-idxNA]
            } else {
                binDatTCNwithoutNA <- binDatTCN
                bkpPosByCHR[[cc]] <- bkpPos
            }
            ## Define Y1 and Y2
            Y <- rbind(Y, binDatTCNwithoutNA)
            colnames(Y) <- names
            Y1 <- NA_real_
            Y2 <- NA_real_
        }
        bkpPosByCHR[[cc]] <- c(min(pos), bkpPosByCHR[[cc]], max(pos))
    }

    structure(list(
      Y1=Y1,
      Y2=Y2,
      Y=Y,
      bkp=bkpPosByCHR
    ), class = c("C3coSegmentation"))
}


#' @export
nbrOfSegments <- function(x, ...) UseMethod("nbrOfSegments")

#' @export
nbrOfSegments.C3coSegmentation <- function(x, ...) {
  nrow(x$Y)
}


#' @export
nbrOfSamples <- function(x, ...) UseMethod("nbrOfSamples")

#' @export
nbrOfSamples.C3coSegmentation <- function(x, ...) {
  ncol(x$Y)
}


#' @export
sampleNames <- function(x) NextMethod("sampleNames")

#' @export
sampleNames.C3coSegmentation <- function(x) {
  colnames(x$Y)
}


#' @export
nbrOfChromosomes <- function(x, ...) UseMethod("nbrOfChromosomes")

#' @export
nbrOfChromosomes.C3coSegmentation <- function(x, ...) {
  length(x$bkp)
}


#' @export
trackNames <- function(x) NextMethod("trackNames")

#' @export
trackNames.C3coSegmentation <- function(x) {
  ## FIXME: Unsafe /HB 2018-03-11
  grep("^Y[0-9]*$", names(x), value = TRUE)
}


#' @export
print.C3coSegmentation <- function(x, ...) {
  s <- sprintf("%s: ", class(x)[1])
  n <- nbrOfSamples(x)
  names <- sampleNames(x)
  if (is.null(names)) {
    names <- seq_len(n)
  } else {
    names <- sQuote(names)
  }
  s <- c(s, sprintf(" Samples: [%d] %s", n, hpaste(names)))
  s <- c(s, sprintf(" Number chromosomes: %d", nbrOfChromosomes(x)))
  s <- c(s, sprintf(" Number segments: %d", nbrOfSegments(x)))
  s <- c(s, sprintf(" Method: jointseg::jointSeg"))
  t <- trackNames(x)
  s <- c(s, sprintf(" Dimensions: [%d] %s", length(t), comma(sQuote(t))))
  s <- paste(s, collapse = "\n")
  cat(s, "\n", sep = "")
}
