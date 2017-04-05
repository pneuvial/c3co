#' Cancer subclone inference
#'
#' @param dat A list of data frames for each patient, each of them of the form
#' \describe{
#'   \item{tcn}{Total copy number}
#'   \item{dh}{Mirrored B allele fraction}
#'   \item{pos}{Position on the genome}
#'   \item{chr}{Chromosome}
#'   }
#' @param parameters.grid A list composed of two vectors named \code{lambda1} and \code{lambda2} of real numbers which are the penalty coefficients for the fused lasso on the minor and major copy number dimension and a vector named \code{nb.arch} of integers which is the number of archetypes in the model
#' @param stat TCN or C1C2
#' @param pathSeg Path to the file that contain segmentation, by default \code{NULL}.
#' @param \dots Further arguments to be passed to \code{\link{fitC3co}}
#' @param verbose A logical value indicating whether to print extra information. Defaults to FALSE
#' @return An object of class [\code{\linkS4class{c3coFit}}]
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
#' l1 <- seq(from=1e-6, to=1e-5, length=3)
#' l2 <- seq(from=1e-6, to=1e-5, length=3)
#' parameters.grid <- list(lambda1=l1, lambda2=l2, nb.arch=2:6)
#' res <- c3co(dat, parameters.grid)
#' resC <- c3co(dat, stat="TCN", parameters.grid)
#'
#' @importFrom methods new
#' @export
c3co <- function(dat, parameters.grid=NULL, stat=c("C1C2", "TCN"), pathSeg=NULL, ..., verbose=FALSE) {
    ## Sanity checks
    stat <- match.arg(stat)
    if (!is.null(dat)) {
        if (!is.list(dat)) {
            stop("Argument 'dat' should be a list ")
        }
        checkCols <- lapply(dat, FUN=function(dd) {
            coln <- colnames(dd)
            ecn <- c("tcn", "dh", "pos", "chr") ## expected
            mm <- match(ecn, coln)
            if (any(is.na(mm))) {
                str <- sprintf("('%s')", paste(ecn, collapse="','"))
                stop("Argument 'dat' should contain columns named", str)
            }
        })
        expectedL <- nrow(dat[[1]]) 
        checklength <- lapply(dat, FUN=function(dd) {
            nrowDD <- nrow(dd)
            if (nrowDD != expectedL) {
                stop("Each data.frame in 'dat' should have the same size")
            }
        })
    }
    if (stat=="TCN") {
        new.getZ <- FALSE
    }
    if (is.null(parameters.grid)) {
        lambda1 <- seq(from=1e-6, to=1e-4, length.out=10)
        lambda2 <- seq(from=1e-6, to=1e-4, length.out=10)
        nb.arch  <- 2:(length(dat)-1)
        parameters.grid <- list(lambda1=lambda1, lambda2=lambda2, nb.arch=nb.arch)
    }
    
    checkGrid <- lapply(names(parameters.grid), FUN=function(na) {
        ecn <- c("lambda1", "lambda2", "nb.arch") ## expected
        mm <- match(na, ecn)
        if (any(is.na(mm))) {
            str <- sprintf("('%s')", paste(ecn, collapse="','"))
            stop("Argument 'parameters.grid' should contain ", str)
        }
    })
    
    if (!is.null(pathSeg)) {
        if (verbose) {
            print("Reading segmentation results from file: ")
            print(pathSeg)
        }
        seg <- readRDS(pathSeg)
    } else {
        seg <- segmentData(dat, stat=stat, verbose=verbose)
    }
    bkpList <- seg$bkp
    
    reslist <- new("c3coFit")
    reslist@bkp <- bkpList
    reslist@segDat <- list(Y1=seg$Y1, Y2=seg$Y2, Y=seg$Y)
    
    if (stat=="C1C2") {
        Y1 <- t(seg$Y1)
        Y2 <- t(seg$Y2)
    } else if (stat=="TCN") {
        Y1 <- t(seg$Y)
        Y2 <- NULL
    }
    
    reslist@fit <- fitC3co(Y1, Y2=Y2, parameters.grid=parameters.grid, ..., verbose=verbose)
    return(reslist)
}
