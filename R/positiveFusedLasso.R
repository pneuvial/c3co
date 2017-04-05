#' Positive fused lasso function
#'
#' @export
#' @param Y1 A matrix containing the segmented minor copy number (\code{n} patients in row and \code{L} segments in columns)
#' @param Y2 A matrix containing the segmented major copy number (\code{n} patients in row and \code{L} segments in columns)
#' @param Z1 A \code{L} x \code{p} matrix containing the \code{L} minor copy numbers of the \code{p} initial latent feature estimates
#' @param Z2 A \code{L} x \code{p} matrix containing the \code{L} major copy numbers of the \code{p} initial latent feature estimates
#' @param lambda1 A real number, the coefficient for the fused penalty for minor copy numbers 
#' @param lambda2 A real number, the coefficient for the fused penalty for major copy numbers 
#' @param eps criterion to stop algorithm (when W do not change sqrt(sum((W-W.old)^2)<eps) 
#' @param max.iter maximum number of iterations of the algorithm
#' @param warn issue a warning if Z1<=Z2 is not always satisfied? Defaults to FALSE
#' @param verbose if you want to print some information during running
#' @return An object of class [\code{\linkS4class{posFused}}]
#' @examples
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 3
#' bkps <- list(c(100,250)*10, c(150,400)*10,c(150,400)*10)
#' regions <-list(c("(0,3)", "(0,1)","(1,2)"), 
#' c("(1,1)", "(0,1)","(1,1)"), c("(0,2)", "(0,1)","(1,1)"))
#' datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN, nbClones, bkps, regions)
#' M <- rbind(c(40, 30, 0), 
#'   c(0, 70, 15),
#'   c(10, 0, 35),
#'   c(15, 0, 0),
#'   c(0, 0, 0))
#' simu <- mixSubclones(subClones=datSubClone, M)
#' seg <- segmentData(simu)
#' lambda <- 0.001
#' Z <- initializeZ(seg$Y1, seg$Y2, nb.arch=4)
#' res <- positiveFusedLasso(seg$Y1, seg$Y2, Z$Z1, Z$Z2, lambda1=lambda, lambda2=lambda)
#' showPosFused(res)
#' resC <- positiveFusedLasso(seg$Y1+seg$Y2, Y2=NULL, Z1=Z$Z, Z2=NULL, lambda1=lambda, lambda2=lambda)
#' showPosFused(resC)
#' 
#' @importFrom methods new
positiveFusedLasso <- function(Y1, Y2, Z1, Z2, lambda1, lambda2, eps=1e-2, max.iter=50, warn=FALSE, verbose=FALSE) {
    n <- nrow(Y1) # number of individuals
    L <- ncol(Y1) # number of loci/segments
    stopifnot(eps>0)
    nb.arch <- ncol(Z1)
    Z <- Z0 <- list(Z1=Z1, Z2=Z2)
    ## Z0 <- list(Z=Z1+Z2, Z1=Z1, Z2=Z2)
    
    iter <- 0
    cond <- FALSE
    delta <- Inf
    ## __________________________________________________
    ## main loop for alternate optimization
    Ymat <- cbind(Y1, Y2)
    lst <- list(Z1=list(Y=Y1, lambda=lambda1))
    if (!is.null(Y2)) {
        lst[["Z2"]] <- list(Y=Y2, lambda=lambda2)
    }

    while (!cond) {
        iter <- iter + 1
        ## __________________________________________________
        ## STEP 1: optimize wrt W (fixed Z1, Z2)
        Zmat <- rbind(Z$Z1, Z$Z2)
        W <- get.W(Zmat, Ymat)
        
        ## __________________________________________________
        ## STEP 2: optimize wrt Z (fixed W)
        Z <- lapply(lst, FUN=function(ll) {  ## TODO: use future_lapply!
            get.Z(ll[["Y"]], ll[["lambda"]], W=W)
        })
       
        ## __________________________________________________
        ## STEP 3: check for convergence of the weights
        if (iter>1) {
            delta <- sqrt(sum((W-W.old)^2))
        }
        cond <- (iter>max.iter || delta<eps)
        
        if (verbose) message("delta:", round(delta, digits=4))
        W.old <- W
    }
    if (verbose) message("Stopped after ", iter, " iterations")
    if (verbose) message("delta:", round(delta, digits=4))
    ## reshape output
    Z1 <- Z$Z1
    Z2 <- Z$Z2
    if (!is.null(Y2)){
        Y <- Y1 + Y2
        Z <- Z1 + Z2
        Y1hat <- W %*% t(Z1)
        Y2hat <- W %*% t(Z2)
    } else {
        Y <- Y1
        Z <- Z1
        Y1hat <- W %*% t(Z)
        Y2hat <- NULL
    }
    
    ## Calculate model fit statistics
    resVar <- sum((Y-W %*% t(Z))^2)
    predVar <- sum((Y-rowMeans(Y))^2)
    loss <- resVar/(n*L)
    kZ <- sum(apply(Z, MARGIN=2L, FUN=diff) != 0)
    BIC <-  n*L*log(loss) + kZ*log(n*L)
    PVE <- 1-resVar/predVar
    
    if (!is.null(Y2) & warn) {  ## sanity check: minor CN < major CN
        dZ <- Z2-Z1
        tol <- 1e-2  ## arbitrary tolerance...
        if (min(dZ) < - tol) {
            warning("For model with ", nb.arch, " features, some components in minor latent profiles are larger than matched components in major latent profiles")
        }
    }
    S <- list(Z=Z, Z1=Z1, Z2=Z2)
    E <- list(Y1=Y1hat, Y2=Y2hat)
    param <- list(nb.arch=nb.arch, lambda1=lambda1, lambda2=lambda2)
    objRes <- new("posFused", S=S, S0=Z0, W=W, E=E, BIC=BIC, PVE=PVE, param=param)
    return(objRes)
}
