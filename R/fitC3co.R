
#' c3co estimation from segment-level copy number data
#' 
#' Estimate c3co model parameters from segment-level copy number data
#' 
#' @param Y1,Y2 Numeric n-by-J matrices containing segment-level copy numbers,
#' where n is the number of samples and J the number of segments.  If `Y2`,
#' which is optional, is specific, then (Y1,Y2) corresponds to minor and major
#' copy numbers, otherwise (Y1) corresponds to total copy numbers.
#'   
#' @param parameters.grid A list composed of two vectors named `lambda1`
#'   and `lambda2` of real numbers which are the penalty coefficients for 
#'   the fused lasso on the minor and major copy number dimension and a vector 
#'   named `nb.arch` of integers which is the number of archetypes in the 
#'   model.
#'   
#' @param warn If `TRUE` and `Y1` is specified, then a warning is produced
#'   if \eqn{Z1 <= Z2} is not satisfied for a candidate number of subclones.
#'   
#' @param intercept logical: should an intercept be included in the model?
#' Defaults to `FALSE`.
#' 
#' @param \dots Further arguments to be passed to [initializeZt()].
#'   
#' @param verbose A logical indicating whether to print extra information.
#'   Defaults to `FALSE`.
#'   
#' @return A list of k objects of class [c3coFit][c3coFit-class],
#'   where k is the number of candidate number of subclones.
#'   
#' @examples
#' ## Simulate locus-level (C1,C2) copy-number data
#' set.seed(7)
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10  ## Number of loci
#' K <- 3L        ## Number of subclones
#' n <- 14L       ## Number of samples
#' bkps <- list(c(100, 250)*10, c(150, 400)*10, c(150, 400)*10)
#' regions <- list(c("(0,3)", "(0,2)", "(1,2)"),
#' c("(1,1)", "(0,1)", "(1,1)"), c("(0,2)", "(0,1)", "(1,1)"))
#' datSubClone <- buildSubclones(len, K, bkps, regions, dataAnnotTP, dataAnnotN)
#' W <- rSparseWeightMatrix(nb.samp=n, nb.arch=K, sparse.coeff=0.7)
#' dat <- mixSubclones(subClones=datSubClone, W=W)
#'
#' ## Segment the copy-number data
#' seg <- segmentData(dat)
#'
#' ## Fit C3CO model
#' l1 <- seq(from=1e-8, to=1e-5, length.out=5L)
#' parameters.grid <- list(lambda1=l1, nb.arch=1:3)
#' fitList <- fitC3co(t(seg$Y1), t(seg$Y2), parameters.grid=parameters.grid)
#' ## FIXME (cf Issue #67):
#' ## fitListC <- fitC3co(t(seg$Y), parameters.grid=parameters.grid) 
#' 
#' 
#' ## A simpler example with toy data
#' K <- 3L   ## Number of subclones
#' J <- 6L   ## Number of segments
#' n <- 20L  ## Number of samples
#'
#' l1 <- 1e-4
#' candP <- 2:J
#' parameters.grid <- list(lambda=l1, nb.arch=candP)
#' 
#' dat <- getToyData(n=n, len=100L, nbClones=K, nbSegs=J, eps=0.2)  ## almost noiseless!
#' sdat <- dat$segment
#' 
#' res <- fitC3co(sdat$Y, parameters.grid=parameters.grid)
#' pvePlot2(res$config$best, ylim=c(0.6, 1))
#' mus <- sapply(res$fit, FUN=function(x) x@mu)
#' ## the mus are all identical??
#' 
#' @importFrom matrixStats colMaxs
#' @export
fitC3co <- function(Y1, Y2 = NULL, parameters.grid = NULL, warn = TRUE, intercept = FALSE, ..., verbose = FALSE) {

    ## preparing data
    Y <- list(Y1=Y1)
    if (!is.null(Y2)) Y$Y2 <- Y2
    ## Sanity checks
    stop_if_not(length(unique(lapply(Y, FUN = dim))) == 1L) # are all the dimension equal?

    if (verbose) mprintf("fitC3co() ...\n")
  
    ## problem dimension
    n <- nrow(Y1)   ## Number of samples
    J <- ncol(Y1)   ## Number of segments
    M <- length(Y)  ## Number of signal dimensions
    
    ### Define grids
    parameters <- checkParams(parameters.grid=parameters.grid, M=M, J=J, verbose=verbose)
    configs <- parameters$configs
    Ks <- parameters$nb.arch
    
    ## sanity checks
    if (max(Ks) > J) {
      stop("Cannot fit model(s) where the number of latent features (K) is larger than the number of segments (J=", J, ")")
    }
    
    fitList <- allRes <- allLoss <- vector("list", length = length(Ks))
    bestConfigp <- allConfig <- NULL
    for (ii in seq_along(Ks)) {
        K_ii <- Ks[ii]
        if (verbose) mprintf(" - Iteration #%d (%d latent features) of %d ...\n", ii, K_ii, length(Ks))

        ## Initialization
        bestRes <- NULL
        bestBIC <- -Inf
        bestConfig <- aConf <- NULL
        ## Z0 are now initialized on the uncentered version of the data Issue #65 
        Z0t <- initializeZt(Y$Y1, Y$Y2, K=K_ii, ...)
        if (verbose) {
            mprintf("   + Parameter configuration: (%s)\n",
                    comma(colnames(configs)))
        }
        allRes[[K_ii]] <- list()
        allLoss[[K_ii]] <- list()
        for (cc in seq_len(nrow(configs))) {
            cfg <- configs[cc, , drop = TRUE]
            if (verbose) {
              mprintf("   + configuration #%d (%s) of %d: ",
                  cc, 
                  comma(sprintf("%s = %g", names(cfg), cfg)),
                  nrow(configs))
            }

              
            ## THIS is the function that costs comme computation 
            if (verbose) mprintf("fused lasso")
            if (is.null(Y2)) cfg <- cfg["lambda1"]  ## Should this be an error?
            res <- positiveFusedLasso(Y = Y, Zt = Z0t, lambda = cfg, intercept = intercept)
            
            if (verbose) mprintf(", BIC = ")
            stats <- modelFitStatistics(Y = Reduce(`+`, Y), Yhat = res@E$Y,
                                        What = res@W, Zhatt = res@Zt$Z)
            BIC <- stats[["BIC"]]
            if (verbose) mprintf("%g", BIC)
            aConf <- c(K_ii, cfg, stats[["PVE"]], BIC, stats[["logLik"]], stats[["loss"]])
            ## replace above line by: aConf <- stats
            allConfig <- rbind(allConfig, aConf)
            allRes[[K_ii]][[cc]] <- res
            allLoss[[K_ii]][[cc]] <- stats[["loss"]]
            
            ## Handle when BIC is NA /HB 2019-03-02
            if (is.finite(BIC) && BIC > bestBIC) { ## BIC has improved: update best model
                bestRes <- res
                bestBIC <- BIC
                bestConfig <- aConf
                if (verbose) mprintf("*")
            }
            if (verbose) mprintf("\n")
        }

        ## Handle when BIC is NA /HB 2019-03-03
        if (!is.null(bestRes)) {
          fitList[[ii]] <- bestRes
          bestConfigp <- rbind(bestConfigp, bestConfig)
	} else {
          bestConfigp <- rbind(bestConfigp, NA)
	}
	
        ## sanity check: minor CN < major CN in the best parameter
        ## configurations (not for all configs by default)
        if (!is.null(Y2) && warn) {
            Z <- bestRes@Zt
            dZ <- Z$Z2 - Z$Z1
            tol <- 1e-2  ## arbitrary tolerance...
            if (min(dZ) < -tol) {
                warning("For model with ", K_ii, " features, some components in minor latent profiles are larger than matched components in major latent profiles")
            }
        }

        if (verbose) mprintf(" - Iteration #%d (%d latent features) of %d ... DONE\n", ii, K_ii, length(Ks))
    } ## for (ii ...)
    
    names(fitList) <- Ks
    cns <- c("nb.feat", colnames(configs), "PVE", "BIC", "logLik", "loss")
    
    bestConfigp <- as.data.frame(bestConfigp)
    colnames(bestConfigp) <- cns
    rownames(bestConfigp) <- NULL
    
    allConfig <- as.data.frame(allConfig)
    colnames(allConfig) <- cns
    rownames(allConfig) <- NULL
    configList <- list(best=bestConfigp, all=allConfig, res=allRes, loss=allLoss)

    res <- list(fit=fitList, config=configList)
  
    if (verbose) mprintf("fitC3co() ... DONE\n")
    
    res
}


checkParams <- function(parameters.grid, M, J, verbose) {
    lambda <-  parameters.grid$lambda
    lambda1 <- parameters.grid$lambda1
    lambda2 <- parameters.grid$lambda2
    if (!is.null(lambda)) {
      lambda1 <- lambda
      if (M > 1L) {
          lambda2 <- lambda
          if (verbose) {
            message(" - Only one regularization parameter is provided. Using the same value for lambda[1] and lambda[2]: ", comma(lambda))
          }
          configs <- cbind(lambda1 = lambda1, lambda2 = lambda2)
      } else {
        configs <- cbind(lambda1 = lambda1)
      }
    } else {
      ## Case C1-C2
      if (M == 2L) {
        if (is.null(lambda1) && is.null(lambda2)) {
          lambda <- 10^(-seq(from=6, to=4, length.out=10L))
          if (verbose) {
            message("Regularization parameter lambda is not provided. Using default value: ")
            mstr(lambda)
          }
          configs <- cbind(lambda1 = lambda, lambda2 = lambda)
        } else {
          configs <- expand.grid(lambda1 = lambda1, lambda2 = lambda2) 
        }
        
      ## Case TCN
      } else {
        lambda <- c(lambda1, lambda2)
        if (!is.null(lambda)) {
          if (verbose) {
            message("Only regularization parameter lambda[1] or lambda[2] is used: ", comma(lambda))
          }
        } else {
          lambda <- 10^(-seq(from=6, to=4, length.out=10L))
          if (verbose) {
            message("Regularization parameter lambda or lambda[1] or lambda[2] is not provided. Using value: ", comma(lambda))
          }
        }
        configs <- cbind(lambda1 = lambda)
      }
    }
    
    ## candidate number of subclones
    nb.arch <- parameters.grid$nb.arch
    if (is.null(nb.arch)) {
        nb.arch <- seq(from=2L, to=J-1L, by=1L)
        if (verbose) {
            message("Parameter 'nb.arch' not provided. Using default value: ", comma(nb.arch))
        }
    }
  
    list(nb.arch = nb.arch, configs = configs)
}
