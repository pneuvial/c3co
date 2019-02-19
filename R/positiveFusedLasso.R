#' Positive fused Lasso function
#'
#' @param Y A list of one or two matrices named `Y1` and `Y2` containing the
#' segmented minor (and possibly the major) copy numbers 
#' (n patients in rows and L segments in columns).
#'
#' @param Zt A list of one or two L-by-p matrices names `Z1` and `Z2`
#' containing the L minor (and possibly the major) copy numbers of the
#' p initial latent feature estimates.
#' _Warning_: These matrices are actually the transpose of Z1 and Z2,
#' which is a notation mistake that will be fixed in a future release.
#'
#' @param lambda A numeric with one or two real numbers, the coefficients
#' for the fused penalty for minor (and possibly the major) copy numbers.
#'
#' @param eps Criterion to stop algorithm
#' (when W do not change sqrt(sum((W-W.old)^2) < eps).
#'
#' @param max.iter Maximum number of iterations of the algorithm.
#'
#' @param warn Issue a warning if `all(Z1 <= Z2)` is not always satisfied?
#' Defaults to `FALSE`.
#'
#' @param verbose If you want to print some information during running.
#'
#' @return An object of class [posFused][posFused-class].
#'
#' @examples
#' set.seed(7)
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
#' W <- rSparseWeightMatrix(5, nbClones, 0.7)
#' datList <- mixSubclones(subClones=datSubClone, W)
#' 
#' segData.TCN <- segmentData(datList,"TCN")
#' Y1 <- t(segData.TCN$Y)
#' Y <- list(Y1 = Y1)
#' Z0t.TCN <- initializeZt(Y1, p = 2, flavor = "nmf")
#' Zt <- list(Z1 = Z0t.TCN$Z1)
#' posFused <- positiveFusedLasso(Y, Zt, lambda=1e-3, verbose=TRUE)
#' modelFitStats(posFused)
#' 
#' segData.C1C2 <- segmentData(datList,"C1C2")
#' Y1 <- t(segData.C1C2$Y1)
#' Y2 <- t(segData.C1C2$Y2)
#' Y <- list(Y1 = Y1, Y2 = Y2)
#' Z0t.C1C2 <- initializeZt(Y1, Y2, p=2, flavor = "nmf")
#' Zt <- list(Z1 = Z0t.C1C2$Z1, Z2 = Z0t.C1C2$Z2)
#' posFusedC <- positiveFusedLasso(Y, Zt, lambda=c(1e-3, 1e-3), verbose=TRUE)
#' modelFitStats(posFusedC)
#' 
#' @importFrom methods new
#' @export
positiveFusedLasso <- function(Y, Zt, lambda, eps=1e-1,
                               max.iter=50, warn=FALSE, verbose=FALSE) {
  ## Argument 'Y':
  stop_if_not(is.list(Y))
  M <- length(Y)
  stop_if_not(M >= 1L)
  n <- nrow(Y[[1]])   ## number of samples
  L <- ncol(Y[[1]])   ## number of segments
  
  ## Argument 'Zt':
  stop_if_not(is.list(Zt), length(Zt) == M, nrow(Zt[[1]]) == L)
  p <- ncol(Zt[[1]])  ## number of subclones/archetypes/latent features
  if (p > n) {
    warning("Under-identified problem: more latent features than samples")
  }

  if (M >= 2) {
    for (mm in 2:M) {
      stop_if_not(nrow(Y[[mm]]) == n, ncol(Y[[mm]]) == L,
                  ncol(Zt[[mm]]) == p, nrow(Zt[[mm]]) == L)
    }
  }

  ## Argument 'lambda':
  stop_if_not(is.numeric(lambda), length(lambda) == M, !anyNA(lambda),
              all(lambda >= 0))

  ## Argument 'eps':
  stop_if_not(is.numeric(eps), length(eps) == 1L, !is.na(eps), eps > 0)

  ## Argument 'max.iter':
  stop_if_not(is.numeric(max.iter), length(max.iter) == 1L,
              !is.na(max.iter), max.iter > 0L)


  Yc <- do.call(cbind, args = Y) # stacked signals

  ## __________________________________________________
  ## main loop for alternate optimization
  iter <- 0; cond <- FALSE; delta <- Inf
  while (!cond) {
    iter <- iter + 1
    
    ## __________________________________________________
    ## STEP 1: optimize w.r.t. W (fixed Z)
    ## if (rank of W) < ncol(Zt), there are too many archetypes...
    WtWm1 <- NULL
    while (is.null(WtWm1)) {
      
      ## solve in W (here individuals - i.e. rows of Yc - are independent)
      W <- get.W(Zt = do.call(rbind, args = Zt), Y = Yc)

      ## Check rank deficiency
      QR.W <- qr(W)
      if (QR.W$rank < p) {
        message("W is rank deficent. Removing a latent feature")
### JC: this means that the column of one must be the first column
### if other rank defiency occurs, we remove the first one arbitrarily
##  FIXME: /HB 2019-02-19
        Zt <- lapply(Zt, FUN = function(z) z[,-1])
        ## Remove matched W.old column
        W.old <- W[,-1] 
        p <- p-1
      } else {
        ## use QR decomposition to save time inverting WtW
        WtWm1 <- tcrossprod(backsolve(qr.R(QR.W), x = diag(p)))
      }
    }    

    ## __________________________________________________
    ## STEP 2: optimize w.r.t. Z (fixed W)
    Zt <- mapply(FUN = get.Zt, Y = Y, lambda = lambda, MoreArgs = list(W = W, WtWm1 = WtWm1), SIMPLIFY = FALSE)
    
    ## __________________________________________________
    ## STEP 3: check for convergence of the weights
    if (iter > 1) { delta <- sqrt(sum((W - W.old)^2)) }
    cond <- (iter > max.iter || delta < eps)
    if (verbose) message("delta:", round(delta, digits=4))
    W.old <- W
  }
  if (verbose) message("Stopped after ", iter, " iterations")
  if (verbose) message("delta:", round(delta, digits=4))

  if (length(Y) > 1 && warn) { ## sanity check: minor CN < major CN
    dZt <- Reduce(`-`, rev(Zt))
    tol <- 1e-2  ## arbitrary tolerance...
    if (min(dZt) < -tol) {
       warning("For model with ", p, " features, some components in minor latent profiles are larger than matched components in major latent profiles")
    }
  }
  
  ## reshape output
  ## accumulate Y, Yhat and Z to get the sum of the two clones (used by Morgane in her representation)
  names(Y) <- paste0("Y", 1:M)
### JC: having a list whose first element has the same name is rather ugly
### and not helful at all to the user  
  Y$Y <- Reduce(`+`, Y)
### JC: should be a method  
  Yhat <- lapply(Zt, FUN = function(Zt_) W %*% t(Zt_))
  names(Yhat) <- paste0("Y", 1:M)
  Yhat$Y <- Reduce(`+`, Yhat)

  names(Zt) <- paste0("Z", 1:M)
### JC: same remark than for Y$Y
  Zt$Z <- Reduce(`+`, Zt)

### JC: Useless ???
  names(lambda) <- paste0("lambda", 1:M)
  params <- c(nb.feat=p, lambda)

  ## Sanity checks
  # FIXME: M + 1L because also Y = Y1 + Y2
  stop_if_not(is.list(Y), length(Y) == M + 1L)
  stop_if_not(is.list(Yhat), length(Yhat) == M + 1L)
  stop_if_not(is.list(Zt), length(Zt) == M + 1L)
  stop_if_not(is.matrix(W), nrow(W) == n, ncol(W) == p)
  stop_if_not(is.numeric(params))
  for (mm in 1:(M+1L)) {
    stop_if_not(nrow(Y[[mm]]) == n, ncol(Y[[mm]]) == L)
    stop_if_not(nrow(Yhat[[mm]]) == n, ncol(Yhat[[mm]]) == L)
    stop_if_not(ncol(Zt[[mm]]) == p, nrow(Zt[[mm]]) == L)
  }
  
### JC: why Z is called S outside of this function
### why not calling Y Z and W by their true name like, 
### signals, archetypes, weights, when outside of this function?
  new("posFused", Y=Y, S=Zt, W=W, E=Yhat, params=params)
}

