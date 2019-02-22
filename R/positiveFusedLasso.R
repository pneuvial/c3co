#' Positive fused Lasso function
#'
#' @param Y A list of one or two matrices named `Y1` and `Y2` containing the
#' segmented minor (and possibly the major) copy numbers 
#' (n patients in rows and J segments in columns).
#'
#' @param Zt A list of one or two J-by-K matrices names `Z1` and `Z2`
#' containing the J minor (and possibly the major) copy numbers of the
#' K initial latent feature estimates.
#' _Warning_: These matrices are actually the transpose of Z1 and Z2,
#' which is a notation mistake that will be fixed in a future release.
#'
#' @param lambda A numeric with one or two real numbers, the coefficients
#' for the fused penalty for minor (and possibly the major) copy numbers.
#'
#' @param intercept logical: should an intercept be included in the model.
#' Defaults to `FALSE`.
#' 
#' @param eps Criterion to stop algorithm
#' (when W do not change sqrt(sum((W-W.old)^2) <= eps).
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
#'
#' len <- 700*10  ## Number of loci
#' K <- 2L        ## Number of subclones
#' n <- 5L        ## Number of samples
#'
#' bkps <- list(
#'     c(20, 30, 50, 60)*100, 
#'     c(10, 40, 60)*100)
#' regions <- list(
#'     c("(1,1)", "(1,2)", "(1,1)", "(0,1)", "(1,1)"),
#'     c("(0,1)", "(1,1)", "(1,2)", "(1,1)"))
#' 
#' datSubClone <- buildSubclones(len, K, bkps, regions, eps=0.2)
#' W <- rSparseWeightMatrix(nb.samp=n, nb.arch=K, sparse.coeff=0.7)
#' datList <- mixSubclones(subClones=datSubClone, W=W)
#' 
#' segData.TCN <- segmentData(datList, "TCN")
#' Y1 <- t(segData.TCN$Y)
#' Y <- list(Y1 = Y1)
#' Z0t.TCN <- initializeZt(Y1, K=K, flavor="nmf")
#' Zt <- list(Z1 = Z0t.TCN$Z1)
#' posFused <- positiveFusedLasso(Y, Zt, lambda=1e-3, verbose=TRUE)
#' modelFitStats(posFused)
#' 
#' segData.C1C2 <- segmentData(datList,"C1C2")
#' Y1 <- t(segData.C1C2$Y1)
#' Y2 <- t(segData.C1C2$Y2)
#' Y <- list(Y1 = Y1, Y2 = Y2)
#' Z0t.C1C2 <- initializeZt(Y1, Y2, K=K, flavor="nmf")
#' Zt <- list(Z1 = Z0t.C1C2$Z1, Z2 = Z0t.C1C2$Z2)
#' posFusedC <- positiveFusedLasso(Y, Zt, lambda=c(1e-3, 1e-3), verbose=TRUE)
#' print(posFusedC)
#' modelFitStats(posFusedC)
#' 
#' @importFrom methods new
#' @export
positiveFusedLasso <- function(Y, Zt, lambda, intercept = FALSE, eps = 1e-1,
                               max.iter = 50L, warn = FALSE, verbose = FALSE) {
  ## Argument 'Y':
  stop_if_not(is.list(Y))
  M <- length(Y)      ## Number of signal dimensions
  stop_if_not(M >= 1L)
  n <- nrow(Y[[1]])   ## Number of samples
  J <- ncol(Y[[1]])   ## Number of segments
  
  ## Argument 'Zt':
  stop_if_not(is.list(Zt), length(Zt) == M, nrow(Zt[[1]]) == J)
  K <- ncol(Zt[[1]])  ## number of subclones/archetypes/latent features
  
  if (M >= 2L) {
    for (mm in 2:M) {
      stop_if_not(nrow(Y[[mm]]) == n, ncol(Y[[mm]]) == J,
                  ncol(Zt[[mm]]) == K, nrow(Zt[[mm]]) == J)
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


  ## Is the problem identifiable?
  if (K > n) {
    stop(sprintf("Under-identified problem: more latent features (K = %d) than samples (n = %d)", K, n))
  }

  ## handling intercept terme
  if (intercept) {
    Y_bar  <- lapply(Y , rowMeans)
    Zt_bar <- lapply(Zt, colMeans)
    Yc <- mapply(sweep, x = Y , STATS = Y_bar , MoreArgs = list(MARGIN = 1, FUN  = "-"), SIMPLIFY = FALSE)
    Zt <- mapply(sweep, x = Zt, STATS = Zt_bar, MoreArgs = list(MARGIN = 2, FUN  = "-"), SIMPLIFY = FALSE)
  } else {
    Yc <- Y
  }
  
  ## __________________________________________________
  ## main loop for alternate optimization
  iter <- 1L
  converged      <- FALSE
  rank_deficient <- FALSE
  lsei_failure   <- FALSE
  delta <- Inf
  while (!converged && iter <= max.iter && !rank_deficient && !lsei_failure) {
    ## __________________________________________________
    ## STEP 1: optimize w.r.t. W (fixed Z)
    
    ## solve in W (here individuals - i.e. rows of Yc - are independent)
    W <- get.W(Zt = do.call(rbind, args = Zt), Y = do.call(cbind, args = Yc))
    if (anyNA(W)) {
      message("No solution found in constrained least-squared problem.")
      lsei_failure <- TRUE
      break
    } 

    # Check rank deficiency
    QR.W <- qr(W)
    if (QR.W$rank < K) {
      message("W is rank deficient: there are too many archetypes")
      rank_deficient <- TRUE
      W <- matrix(NA_real_, n, K)
      break
    } else {
      ## use QR decomposition to save time inverting WtW
      WtWm1 <- tcrossprod(backsolve(qr.R(QR.W), x = diag(K)))
    }
    
    ## __________________________________________________
    ## STEP 2: optimize w.r.t. Z (fixed W)
    Zt <- mapply(FUN = get.Zt, Y = Yc, lambda = lambda, MoreArgs = list(W = W, WtWm1 = WtWm1), SIMPLIFY = FALSE)
    if (intercept) { # 
      Zt_bar <- lapply(Zt, colMeans)
      Zt <- mapply(sweep, x = Zt, STATS = Zt_bar, MoreArgs = list(MARGIN = 2, FUN  ="-"), SIMPLIFY = FALSE)
    }
    
    ## __________________________________________________
    ## STEP 3: check for convergence of the weights
    if (iter > 1L) {
      delta <- sqrt(sum((W - W.old)^2))
      converged <- (delta <= eps)
    }
    if (verbose) message("delta:", round(delta, digits=4L))

    ## Next iteration
    W.old <- W
    iter <- iter + 1L
  } ## while (...)
  
  if (verbose) {
    if (converged) {
      message("Converged after ", iter, " iterations")
    } else if (rank_deficient) {
      message("Stopped after ", iter, " because of rank deficiency")
    } else if (lsei_failure) {
      message("Stopped after ", iter, " because of lsei fails")
    } else {
      message("Stopped after ", iter, " iterations without reaching convergence")
    }
    message("delta:", round(delta, digits=4L))
  }

  if (warn && length(Y) > 1L) { ## sanity check: minor CN < major CN
    dZt <- Reduce(`-`, rev(Zt))
    tol <- 1e-2  ## arbitrary tolerance...
    if (min(dZt) < -tol) {
       warning("For model with ", K, " features, some components in minor latent profiles are larger than matched components in major latent profiles")
    }
  }

  ## list of Intercept terms
  if (intercept) { 
    mu <- mapply(FUN = function(y_bar, zt_bar) {
        as.numeric(y_bar - W %*% zt_bar)
      }, y_bar = Y_bar, zt_bar = Zt_bar, SIMPLIFY = FALSE
    )
  } else {
    mu <- rep(list(rep(0,n)), M)
  }
  
  ## reshape output
  
  ## accumulate Y, Yhat and Z to get the sum of the two clones (used by Morgane in her representation)
  ## JC 2019/02/22: check when this representation is used !!
  
  names(Y) <- paste0("Y", 1:M)
  Y$Y <- Reduce(`+`, Y)
  
  Yhat <- mapply(function(mu_, Zt_) {mu_ + W %*% t(Zt_)}, mu_ = mu, Zt_ = Zt, SIMPLIFY = FALSE)
  names(Yhat) <- paste0("Y", 1:M)
  Yhat$Y <- Reduce(`+`, Yhat)

  names(Zt) <- paste0("Z", 1:M)
  Zt$Z <- Reduce(`+`, Zt)

  ##JC 2019/02/22: do not knwo if accumulated mu is used (should be) when computing repsentation and stats
  names(mu) <- paste0("mu", 1:M)
  mu$mu <- Reduce(`+`, mu)

  names(lambda) <- paste0("lambda", 1:M)
  params <- c(nb.feat = K, lambda)

  ## Sanity checks
  # FIXME: M + 1L because also Y = Y1 + Y2
  stop_if_not(is.list(Y), length(Y) == M + 1L)
  stop_if_not(is.list(mu), length(mu) == M + 1L)
  stop_if_not(is.list(Yhat), length(Yhat) == M + 1L)
  stop_if_not(is.list(Zt), length(Zt) == M + 1L)
  stop_if_not(is.matrix(W), nrow(W) == n, ncol(W) == K)
  stop_if_not(is.numeric(params))
  for (mm in 1:(M + 1L)) {
    stop_if_not(length(mu[[mm]]) == n)
    stop_if_not(nrow(Y[[mm]]) == n, ncol(Y[[mm]]) == J)
    stop_if_not(nrow(Yhat[[mm]]) == n, ncol(Yhat[[mm]]) == J)
    stop_if_not(ncol(Zt[[mm]]) == K, nrow(Zt[[mm]]) == J)
  }
  
#  ## Sanity checks of estimated values (Issue #49)
#  tol <- sqrt(.Machine$double.eps)
#  stop_if_not(all(is.finite(W)), all(W >= 0 - tol), all(W <= 1 + tol))
#  for (mm in 1:(M+1L)) {
#    stop_if_not(
#      all(is.finite(Zt[[mm]])),
#      all(Zt[[mm]] >= 0 - tol)   ## Does not always look to be true
#    )
#  }
  
  new("posFused",
      Y          = Y        ,
      mu         = mu       ,
      Zt         = Zt       ,
      W          = W        ,
      E          = Yhat     ,
      params     = params   ,
      converged  = converged,
      iterations = iter
    )
}
