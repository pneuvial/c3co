#' Cancer subclone inference
#'
#' @param dat A list of data frames for each patient, each of them of the form:
#' \describe{
#'   \item{`tcn`}{Total copy number}
#'   \item{`dh`}{Mirrored B allele fraction}
#'   \item{`pos`}{Position on the genome}
#'   \item{`chr`}{Chromosome}
#' }
#'
## FIXME: What does ".. list composed of : the ..." mean? /HB 2018-03-04
#' @param parameters.grid Is a list composed of : the penalty coefficients
#' named either `lambda`, `lambda1` and `lambda2`. For the model on minor and
#' major CN, it is possible to give two grids (one for minor CN and one for
#' major CN), and all combination is tested. If only one grid is given the
#' consider  that `lambda1` = `lambda2`.
#' and composed of : a vector named `nb.arch` of integers which is the number
#' of features in the model.
#'
#' @param stat `"TCN"` or `"C1C2"`.
#'
#' @param segDat Either a path to the file that contains segmentation
#' (\file{*.rds} file), by default `NULL`, or a file that contains segmentation 
#' (e.g. from [segmentData()]).
#'
#' @param \dots Further arguments to be passed to [fitC3co()].
#'
#' @param verbose A logical value indicating whether to print extra
#' information. Defaults to `FALSE`.
#'
#' @return An object of class [c3coFit][c3coFit-class].
#'
#' @examples
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 3
#' set.seed(88)
#' bkps <- list(c(100, 250)*10, c(150, 400)*10, c(150, 400)*10)
#' regions <- list(c("(1,2)", "(0,2)", "(1,2)"),
#' c("(0,3)", "(0,1)", "(1,1)"), c("(0,2)", "(0,1)", "(1,1)"))
#' datSubClone <- buildSubclones(len, nbClones, bkps, regions, dataAnnotTP, dataAnnotN)
#' M <- rSparseWeightMatrix(15, 3, 0.7)
#' dat <- mixSubclones(subClones=datSubClone, M)
#' l1 <- seq(from=1e-6, to=1e-4, length.out=10)
#' parameters.grid <- list(lambda=l1, nb.arch=2:6)
#' res <- c3co(dat, parameters.grid)
#' l2 <- seq(from=1e-6, to=1e-5, length.out=10)
#' parameters.grid.2 <- list(lambda=l2, nb.arch=2:6)
#' resC <- c3co(dat, stat="TCN", parameters.grid.2)
#' @importFrom methods new
#' @export
c3co <- function(dat, parameters.grid=NULL, stat=c("C1C2", "TCN"),
                 segDat=NULL, ..., verbose=FALSE) {
    ## Sanity checks
    stat <- match.arg(stat)
    if (!is.null(dat)) {
        if (!is.list(dat)) {
            stop("Argument 'dat' should be a list")
        }
        checkCols <- lapply(dat, FUN=function(dd) {
            coln <- colnames(dd)
            ecn <- c("tcn", "dh", "pos", "chr") ## expected
            mm <- match(ecn, coln)
            if (anyNA(mm)) {
                stop("Argument 'dat' should contain columns named ",
                     comma(sQuote(ecn)))
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
    if (is.null(parameters.grid)) {
        lambda1 <- seq(from=1e-6, to=1e-4, length.out=10)
        nb.arch <- 2:(length(dat)-1)
        parameters.grid <- list(lambda=lambda1, nb.arch=nb.arch)
    }

    checkGrid <- lapply(names(parameters.grid), FUN=function(na) {
        ecn <- c("lambda", "lambda1", "lambda2", "nb.arch") ## expected
        mm <- match(na, ecn)
        if (anyNA(mm)) {
            stop("Argument 'parameters.grid' should contain ",
                 comma(sQuote(ecn)))
        }
    })

    if (!is.null(segDat)) {
      if (is.character(segDat)) {
        if (verbose) {
            message("Reading segmentation results from file: ")
            mprint(segDat)
        }
        seg <- readRDS(segDat)
      }else{
        ## Sanity check
        checkGrid <- lapply(names(segDat), FUN=function(na) {
          ecn <- c("bkp", "Y1", "Y2", "Y") ## expected
          mm <- match(na, ecn)
          if (anyNA(mm)) {
              stop("Argument 'parameters.grid' should contain ",
                   comma(sQuote(ecn)))
          }
        })
        message("Segmented data is provided, skip segment step")
        seg <- segDat
      }
    } else {
        seg <- segmentData(dat, stat=stat, verbose=verbose)
    }
    bkpList <- seg$bkp

    reslist <- new("c3coFit")
    reslist@bkp <- bkpList
    reslist@segDat <- list(Y1=seg$Y1, Y2=seg$Y2, Y=seg$Y)

    if (stat == "C1C2") {
        Y1 <- t(seg$Y1)
        Y2 <- t(seg$Y2)
    } else if (stat == "TCN") {
        Y1 <- t(seg$Y)
        Y2 <- NULL
    }

    fit <- fitC3co(Y1, Y2=Y2, parameters.grid=parameters.grid,
                           ..., verbose=verbose)
    reslist@fit <- fit$fit
    reslist@config <- fit$config
    
    reslist
}
