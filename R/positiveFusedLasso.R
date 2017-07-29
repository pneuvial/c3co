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
#' len <- 700*10
#' nbClones <- 2
#' bkps <- list(
#'     c(20, 30, 50, 60)*100, 
#'     c(10, 40, 60)*100)
#' regions <- list(
#'     c("(1,1)", "(1,2)", "(1,1)", "(0,1)", "(1,1)"),
#'     c("(0,1)", "(1,1)", "(1,2)", "(1,1)"))
#' 
#' datSubClone <- buildSubclones(len, nbClones, bkps, regions, eps=0.2)
#' W <-  rSparseWeightMatrix(5, nbClones, 0.7)
#' datList <- mixSubclones(subClones=datSubClone, W)
#' 
#' segData.TCN  <- segmentData(datList,"TCN")
#' Y1 <- t(segData.TCN$Y)
#' Y <- list(Y1 = Y1)
#' Z0.TCN <- initializeZ(Y1, p = 2, flavor = "nmf")
#' Z <- list(Z1 = Z0.TCN$Z1)
#' posFused <- positiveFusedLasso(Y, Z, 1e-3, verbose=TRUE)
#' posFused 
#' 
#' segData.C1C2 <- segmentData(datList,"C1C2")
#' Y1 <- t(segData.C1C2$Y1)
#' Y2 <- t(segData.C1C2$Y2)
#' Y <- list(Y1 = Y1, Y2 = Y2)
#' Z0.C1C2 <- initializeZ(Y1, Y2, p=2, flavor = "nmf")
#' Z <- list(Z1 = Z0.C1C2$Z1,Z2 = Z0.C1C2$Z2)
#' posFusedC <- positiveFusedLasso(Y, Z, c(1e-3, 1e-3), verbose=TRUE)
#' posFusedC
#' 
#' @importFrom methods new
#' @export
positiveFusedLasso <- function(Y, Z, lambda, eps=1e-1,
                               max.iter=50, warn=FALSE, verbose=FALSE) {
  ## problem dimensions
  M <- length(Y)  # number of signal (one or two)
  stopifnot(length(lambda) == M)
  
  Z0 <- Z ## save original Z
  ## Yc is a list of matrices centered row-wise
  Yc   <- lapply(Y, function(y) sweep(y, 1, rowMeans(y), "-"))
  # the vector of means averaged over the M signal (required for the intercept)
  Ybar <- Reduce("+",lapply(Y, rowMeans))/M # average over the signals
  
  ## __________________________________________________
  ## main loop for alternate optimization
  iter <- 0
  cond <- FALSE
  delta <- Inf
  failure <- FALSE
  while (!cond) {
    iter <- iter + 1
    ## __________________________________________________
    ## STEP 1: optimize w.r.t. W (fixed Z)
    
    ## matrices Z must be centered column-wise
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
    if (any(failure)) {
      Z <- Z0
    }
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

  if(M > 1 & warn) { ## sanity check: minor CN < major CN
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
  
  ## reshape output
  Yhat <- lapply(Z, function(Z_) sweep(W %*% t(Z_), 1, mu, "+"))  
  names(Yhat) <- paste0("Y", 1:M)
  ## accumulate Yhat and Z to get the sum of the two clones (used by Morganne in its representation)
  Yhat$Y <- Reduce("+",Yhat)
  names(Z)    <- paste0("Z", 1:M)
  Z$Z    <- Reduce("+",Z)
  
  objRes <- new("posFused", S=Z, S0=Z0, W=W, mu=mu, E=Yhat, failure=failure)
  return(objRes)
}

