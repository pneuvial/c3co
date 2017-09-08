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
#' modelFitStats(posFused)
#' 
#' segData.C1C2 <- segmentData(datList,"C1C2")
#' Y1 <- t(segData.C1C2$Y1)
#' Y2 <- t(segData.C1C2$Y2)
#' Y <- list(Y1 = Y1, Y2 = Y2)
#' Z0.C1C2 <- initializeZ(Y1, Y2, p=2, flavor = "nmf")
#' Z <- list(Z1 = Z0.C1C2$Z1,Z2 = Z0.C1C2$Z2)
#' posFusedC <- positiveFusedLasso(Y, Z, c(1e-3, 1e-3), verbose=TRUE)
#' modelFitStats(posFusedC)
#' 
#' @importFrom methods new
#' @export
positiveFusedLasso <- function(Y, Z, lambda, eps=1e-1,
                               max.iter=50, warn=FALSE, verbose=FALSE) {

  ## problem dimensions
  M <- length(Y)
  p <- nrow(Z[[1]])  ## number of subclones
  S <- ncol(Z[[1]])  ## number of archetypes

### JC: in my opinion, these checks are useless 
### because this function is not meant to be used by the user, right?
  stopifnot(length(lambda) == M)
  stopifnot(length(Z)      == M)

  Yc <- do.call(cbind, Y) # stacked signals

  ## __________________________________________________
  ## main loop for alternate optimization
  iter <- 0; cond <- FALSE; delta <- Inf
  while (!cond) {
    iter <- iter + 1
    
    ## __________________________________________________
    ## STEP 1: optimize w.r.t. W (fixed Z)
    ## if (rank of W) < ncol(Z), there are too many archetypes...
    WtWm1 <- NULL
    while (is.null(WtWm1)) {
      
      ## solve in W (here individuals - i.e. rows of Yc - are independent)
      W <- get.W(do.call(rbind, Z), Yc)

      ## Check rank deficiency
      QR.W <- qr(W)
      if (QR.W$rank < S) {
        message("W is rank deficent. Removing an archetype")
### JC: this means that the column of one must be the first column
### if other rank defiency occurs, we remove the first one arbitrarly
        Z <- lapply(Z, function(z) z[,-1])
        S <-  S-1
      } else {
        ## use QR decomposition to save time inverting WtW
        WtWm1  <- tcrossprod(backsolve(qr.R(QR.W),diag(S)))
      }
    }    

    ## __________________________________________________
    ## STEP 2: optimize w.r.t. Z (fixed W)
    Z <- mapply(get.Z, Y, lambda, MoreArgs = list(W, WtWm1), SIMPLIFY=FALSE)
    
    ## __________________________________________________
    ## STEP 3: check for convergence of the weights
    if (iter > 1) {delta <- sqrt(sum((W - W.old)^2))}
    cond <- (iter > max.iter || delta < eps)
    if (verbose) message("delta:", round(delta, digits=4))
    W.old <- W
  }
  if (verbose) message("Stopped after ", iter, " iterations")
  if (verbose) message("delta:", round(delta, digits=4))

  if(length(Y) > 1 & warn) { ## sanity check: minor CN < major CN
    dZ <- Reduce("-", rev(Z))
    tol <- 1e-2  ## arbitrary tolerance...
    if (min(dZ) < -tol) {
       warning("For model with ", p, " features, some components in minor latent profiles are larger than matched components in major latent profiles")
    }
  }
  
  ## reshape output
  ## accumulate Y, Yhat and Z to get the sum of the two clones (used by Morgane in her representation)
  names(Y) <- paste0("Y", 1:M)
### JC: having a list whose first element has the same name is rather ugly
### and not helful at all to the user  
  Y$Y <- Reduce("+",Y)
### JC: should be a method  
  Yhat <- lapply(Z, function(Z_) W %*% t(Z_))
  names(Yhat) <- paste0("Y", 1:M)
  Yhat$Y <- Reduce("+",Yhat)

  names(Z)    <- paste0("Z", 1:M)
### JC: same remark than for Y$Y
  Z$Z    <- Reduce("+",Z)

### JC: Useless ???
  names(lambda) <- paste0("lambda", 1:M)
  params <- c(nb.feat=p, lambda)

### JC: why Z is called S outside of this function
### why not calling Y Z and W by their true name like, 
### signals, archetypes, weights, when outside of this function?
  objRes <- new("posFused", Y=Y, S=Z, W=W, E=Yhat, params=params)
  return(objRes)
}

