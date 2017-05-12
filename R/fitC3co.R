#' c3co estimation from segment-level copy number data
#' 
#' Estimate c3co model parameters from segment-level copy number data
#' 
#' 
#' @param Y1 A numeric n x S matrix, segment-level minor copy numbers. n is the
#'   number of samples and S the number of segments
#'   
#' @param Y2 An optional numeric n x S matrix, segment-level major copy numbers.
#'   If \code{NULL}, the model is estimated on Y1 only
#'   
#' @param parameters.grid A list composed of two vectors named \code{lambda1} 
#'   and \code{lambda2} of real numbers which are the penalty coefficients for 
#'   the fused lasso on the minor and major copy number dimension and a vector 
#'   named \code{nb.arch} of integers which is the number of archetypes in the 
#'   model
#'   
#' @param warn issue a warning if Z1 <= Z2 is not satisfied for a candidate 
#'   number of subclones? Defaults to TRUE
#'   
#' @param \dots Further arguments to be passed to 
#'   \code{\link{positiveFusedLasso}}
#'   
#' @param verbose A logical value indicating whether to print extra information.
#'   Defaults to FALSE
#'   
#' @return A list of \code{k} objects of class [\code{\linkS4class{c3coFit}}],
#'   where \code{k} is the number of candidate number of subclones
#'   
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
#' M <- rSparseWeightMatrix(15, 3, sparse.coeff=0.7)
#' dat <- mixSubclones(subClones=datSubClone, M)
#' seg <- segmentData(dat)
#' 
#' l1 <- seq(from=1e-6, to=1e-4, length.out=5)
#' parameters.grid <- list(lambda=l1, nb.arch=2:6)
#' fitList <- fitC3co(t(seg$Y1), t(seg$Y2), parameters.grid=parameters.grid, verbose=TRUE)
#' fitListC <- fitC3co(t(seg$Y), parameters.grid=parameters.grid)
#' 
#' @importFrom methods slot
#' @importFrom matrixStats colMaxs
#' @export
fitC3co <- function(Y1, Y2=NULL, parameters.grid=NULL, warn=TRUE,
                    ..., verbose=FALSE) {
    ## Sanity checks
    
    n <- nrow(Y1)
    nseg <- ncol(Y1)
    if (!is.null(Y2)) {
        stopifnot(nrow(Y2) == n)
        stopifnot(ncol(Y2) == nseg)
    }
    lambda <- 10^(-seq(from=6, to=4, length.out=10))
    lambda1 <- parameters.grid$lambda
    if (is.null(lambda1)) {
        lambda1 <- lambda
        if (verbose) {
            message("Regularization parameter lambda[1] not provided. Using default value: ")
            mstr(lambda1)
        }
    }
    
    lambda2 <- parameters.grid$lambda
    configs <- expand.grid(lambda1=lambda1)  ## for when Y2 is NULL
    if (!is.null(Y2)) {
        if (is.null(lambda2)) {
            lambda2 <- lambda
            if (verbose) {
                message("Regularization parameter lambda[2] not provided. Using default value: ")
                mstr(lambda2)
            }
        }
        configs <- cbind(lambda1=lambda1, lambda2=lambda2)
    }
    
    ## candidate number of subclones
    nb.arch <- parameters.grid$nb.arch
    if (is.null(nb.arch)) {
        nb.arch  <- seq(from=2, to=min(nseg, n)-1, by=1)
        if (verbose) {
            message("Parameter 'nb.arch' not provided. Using default value: ")
            mstr(nb.arch)
        }
    }
    rm(list = "parameters.grid")
    
    it <- 1
    pp <- nb.arch[it]
    cond <- FALSE  ## condition for (early) stopping
    fitList <- list()
    while (!cond) {
        if (verbose) {
            message("Number of latent features: ", pp)
        }
        ## Initialization
        BICp <- +Inf
        bestConfig <- NULL
        
        Z0 <- initializeZ(Y1, Y2=Y2, nb.arch=pp, ...)
        
        if (verbose) {
            message("Parameter configuration: (",
                    paste(colnames(configs), collapse="; "), ")")
        }
        for (cc in seq_len(nrow(configs))) {
            cfg <- configs[cc, , drop=FALSE]
            if (verbose) {
                message(paste(cfg, collapse="; "))
            }
            l1 <- cfg[, "lambda1"]
            l2 <- NULL
            if (!is.null(Y2)) l2 <- cfg[, "lambda2"]
            res <- positiveFusedLasso(Y1, Y2=Y2, Z1=Z0$Z1, Z2=Z0$Z2,
                                      lambda1=l1, lambda2=l2, verbose=FALSE)
            if (res@BIC < BICp) { ## BIC has improved: update best model
                res.l <- res
                BICp <- res@BIC
                bestConfig <- cfg
            }
        }
        fitList[[it]] <- res.l
        ## sanity check: minor CN < major CN in the best parameter
        ## configurations (not for all configs by default)
        if (!is.null(Y2) & warn) {
            Z <- slot(res.l, "S")
            dZ <- Z$Z2 - Z$Z1
            tol <- 1e-2  ## arbitrary tolerance...
            if (min(dZ) < - tol) {
                warning("For model with ", pp, " features, some components in minor latent profiles are larger than matched components in major latent profiles")
            }
        }
        
        it <- it + 1
        pp <- nb.arch[it]
        ## stop if pp has reached the max of its grid
        cond <-  is.na(pp)
    }
    return(fitList)
}
