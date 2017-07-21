#' Positive fused lasso function
#'
#' @param Y1 A matrix containing the segmented minor copy number
#' (\code{n} patients in row and \code{L} segments in columns)
#'
#' @param Y2 A matrix containing the segmented major copy number
#' (\code{n} patients in row and \code{L} segments in columns)
#'
#' @param Z1 A \code{L} x \code{p} matrix containing the \code{L}
#' minor copy numbers of the \code{p} initial latent feature estimates
#'
#' @param Z2 A \code{L} x \code{p} matrix containing the \code{L}
#' major copy numbers of the \code{p} initial latent feature estimates
#'
#' @param lambda1 A real number, the coefficient for the fused penalty
#' for minor copy numbers
#'
#' @param lambda2 A real number, the coefficient for the fused penalty
#' for major copy numbers
#'
#' @param eps criterion to stop algorithm (when W do not change
#' sqrt(sum((W-W.old)^2) < eps)
#'
#' @param max.iter maximum number of iterations of the algorithm
#'
#' @param warn issue a warning if Z1 <= Z2 is not always satisfied?
#' Defaults to FALSE
#'
#' @param verbose if you want to print some information during running
#'
#' @return An object of class [\code{\linkS4class{posFused}}]
#'
#' @examples
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 3
#' set.seed(2)
#' bkps <- list(c(100, 250)*10, c(150, 400)*10, c(150, 400)*10)
#' regions <- list(c("(0,3)", "(0,1)", "(1,2)"),
#' c("(1,1)", "(0,1)", "(1,1)"), c("(0,2)", "(0,1)", "(1,1)"))
#' datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN,
#'                               nbClones, bkps, regions)
#' M <- rSparseWeightMatrix(15, 3, 0.4)
#' simu <- mixSubclones(subClones=datSubClone, M)
#' seg <- segmentData(simu)
#' lambda <- 0.01
#' Z <- initializeZ(t(seg$Y1), t(seg$Y2), nb.arch=3)
#' res <- positiveFusedLasso(t(seg$Y1), t(seg$Y2), Z$Z1, Z$Z2,
#'                           lambda1=lambda, lambda2=lambda, verbose=TRUE)
#' res
#' resC <- positiveFusedLasso(t(seg$Y1+seg$Y2), Y2=NULL, Z1=Z$Z, Z2=NULL,
#'                            lambda1=lambda, lambda2=lambda)
#' resC
#'
#' @importFrom methods new
#' @export
positiveFusedLasso <- function(Y1, Y2, Z1, Z2, lambda1, lambda2, eps=1e-1,
                               max.iter=50, warn=FALSE, verbose=FALSE) {
    n <- nrow(Y1) # number of individuals
    L <- ncol(Y1) # number of loci/segments
    stopifnot(eps > 0)
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
        if(iter == 1){
          delta <- tryCatch({
            i <- solve(t(W) %*% W)
            delta <- Inf
          }, error = function(err) {
            # error handler picks up where error was generated
            if(verbose) message("ERROR to solve inv(WtW):  ",err)
            delta <- 0
            return(delta)
          })
          if(delta == 0){
            if(verbose) message("No solution of this combination of lambda")
          }else{
            ## STEP 2: optimize wrt Z (fixed W)
            Z <- lapply(lst, FUN = function(ll) {
              get.Z(Y = ll[["Y"]], lambda = ll[["lambda"]], W=W)
            })
          }
        }else{
          W <- tryCatch({
            i <- solve(t(W) %*% W)
            W
          }, error = function(err) {
            # error handler picks up where error was generated
            if(verbose) message("ERROR to solve inv(WtW):  ",err)
            return(W.old)
          })
          ## STEP 2: optimize wrt Z (fixed W)
          Z <- lapply(lst, FUN = function(ll) {
            get.Z(Y = ll[["Y"]], lambda = ll[["lambda"]], W=W)
          })
          delta <- sqrt(sum((W - W.old)^2))
        }
        
        ## __________________________________________________
        ## STEP 3: check for convergence of the weights
        cond <- (iter > max.iter || delta < eps)
        if (verbose) message("delta:", round(delta, digits=4))
        W.old <- W
    }
    if (verbose) message("Stopped after ", iter, " iterations")
    if (verbose) message("delta:", round(delta, digits=4))
    
    ## reshape output
    Z1 <- Z$Z1
    Z2 <- Z$Z2
    if (!is.null(Y2)) {
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
    if (!is.null(Y2) & warn) {  ## sanity check: minor CN < major CN
        dZ <- Z2-Z1
        tol <- 1e-2  ## arbitrary tolerance...
        if (min(dZ) < - tol) {
            warning("For model with ", nb.arch, " features, some components in minor latent profiles are larger than matched components in major latent profiles")
          idx <- 1:ncol(Z2)
          Z1 <- sapply(idx, function(ii){
            jj <- which(Z1[,ii]>Z2[,ii])
            Z1[jj, ii] <- Z2[jj, ii]
            return(Z1 [,ii])
          })
        }
    }
    S <- list(Z=Z, Z1=Z1, Z2=Z2)
    E <- list(Y1=Y1hat, Y2=Y2hat)
    objRes <- new("posFused",
                  S=S, S0=Z0, W=W, E=E)
    return(objRes)
}
