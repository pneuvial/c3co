#' Positive fused lasso function
#'
#' @param Y A list of one or two matrices named Y1 and Y2 containing the segmented minor 
#' (and possibly the major) copy number 
#' (\code{n} patients in row and \code{L} segments in columns)
#'
#' @param Z A list of one or two \code{L} x \code{p} matrices names Z1 and Z2 containing the \code{L}
#' minor (and possibly the major) copy numbers of the \code{p} initial 
#' latent feature estimates
#'
#' @param lambda A numeric with one or two real numbers, the coefficients for the fused penalty
#' for minor (and possibly the major) copy numbers
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
#' datSubClone <- buildSubclones(len, nbClones, bkps, regions, dataAnnotTP, dataAnnotN)
#' M <- rSparseWeightMatrix(15, 3, 0.4)
#' simu <- mixSubclones(subClones=datSubClone, M)
#' seg <- segmentData(simu)
#' lambda <- 0.01
#' Z <- initializeZ(t(seg$Y1), t(seg$Y2), p=3)
#' res <- positiveFusedLasso(t(seg$Y1), t(seg$Y2), Z$Z1, Z$Z2,
#'                           lambda1=lambda, lambda2=lambda, verbose=TRUE)
#' res
#' resC <- positiveFusedLasso(t(seg$Y1+seg$Y2), Y2=NULL, Z1=Z$Z, Z2=NULL,
#'                            lambda1=lambda, lambda2=lambda)
#' resC
#'
#' @importFrom methods new
#' @export
positiveFusedLasso <- function(Y, Z, lambda, eps=1e-1,
                               max.iter=50, warn=FALSE, verbose=FALSE) {

  ## problem dimensions
  M <- length(Y)  # number of signal (one or two)
  stopifnot(length(lambda) == M)
  Z0 <- Z
  ## Yc is a list of matrices centered row-wise
  Yc   <- lapply(Y, function(y) sweep(y, 1, rowMeans(y), "-"))
  # the vector of means averaged over the M signal (required for the intercept)
  Ybar <- Reduce("+",lapply(Y, rowMeans))/M # average over the signals
  
  ## __________________________________________________
  ## main loop for alternate optimization
  iter <- 0
  cond <- FALSE
  delta <- Inf
  while (!cond) {
    iter <- iter + 1
    ## __________________________________________________
    ## STEP 1: optimize w.r.t. W (fixed Z)
    
    ## Matrices Z must be centering columnwise
    ## the list of achetype matrices centered accordingly
    Zc <- lapply(Z, function(z) sweep(z, 2, colMeans(z), "-"))
    # the vector of means averaged over the M signal (required for the intercept)
    Zbar <- Reduce("+",lapply(Z, colMeans))/M # average of the M signal

    ## solve in W (here individuals - i.e. rows of Yc - are independent)
    W <- get.W(do.call(rbind, Zc), do.call(cbind, Yc))
    
    ## the vector of intercept (one per patient)
    mu <- Ybar - as.numeric(tcrossprod(Zbar, W))

    ## __________________________________________________
    ## STEP 2: optimize w.r.t. Z (fixed W)
    ## centering Y with the current intercept mu
    Yc.mu <- lapply(Y, function(y) sweep(y, 1, mu, "-"))
    Z <- lapply(1:M, function(m) {
      return(get.Z(Yc.mu[[m]], W, lambda[m]))
    })
    failure <- sapply(Z, inherits, "try-error")
    if (any(failure))
      Z <- Z0
    
    if (iter > 1)
      delta <- sqrt(sum((W - W.old)^2))

    ## __________________________________________________
    ## STEP 3: check for convergence of the weights
    cond <- (iter > max.iter || delta < eps || any(failure))
    if (verbose) message("delta:", round(delta, digits=4))
    W.old <- W
  }
  if (verbose) message("Stopped after ", iter, " iterations")
  if (verbose) message("delta:", round(delta, digits=4))

  ## reshape output
  Yhat <- lapply(Z, function(Z_) sweep(W %*% t(Z_), 1, mu, "+"))  
  names(Yhat) <- paste0("Yhat", 1:M)
  names(Z)    <- paste0("Z", 1:M)

  if(length(Yhat) > 1 & warn) { ## sanity check: minor CN < major CN
    dZ <- Reduce("-", rev(Z))
    tol <- 1e-2  ## arbitrary tolerance...
    if (min(dZ) < -tol) {
       warning("For model with ", nb.arch, " features, some components in minor latent profiles are larger than matched components in major latent profiles")
  #     idx <- 1:ncol(Z2)
  #     Z1 <- sapply(idx, function(ii){
  #       jj <- which(Z1[,ii]>Z2[,ii])
  #       Z1[jj, ii] <- Z2[jj, ii]
  #       return(Z1 [,ii])
  #     })
    }
    
  }
  
  ## add Z, the sum of the two clones (use by Morganne in its representation)
  Z$Z <- Reduce("+",Z)
  objRes <- new("posFused", S=Z, S0=Z0, W=W, mu=mu, E=Yhat)
  return(objRes)
}

