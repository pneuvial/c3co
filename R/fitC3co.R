#' c3co estimation from segment-level copy number data
#' 
#' Estimate c3co model parameters from segment-level copy number data
#' 
#' @param Y1 A numeric n-by-S matrix, segment-level minor copy numbers, where n
#'   is the number of samples and S the number of segments.
#'   
#' @param Y2 An optional numeric n-by-S matrix, segment-level major copy numbers.
#'   If `NULL`, the model is estimated on `Y1` only.
#'   
#' @param parameters.grid A list composed of two vectors named `lambda1`
#'   and `lambda2` of real numbers which are the penalty coefficients for 
#'   the fused lasso on the minor and major copy number dimension and a vector 
#'   named `nb.arch` of integers which is the number of archetypes in the 
#'   model.
#'   
#' @param warn Issue a warning if \eqn{Z1 <= Z2} is not satisfied for a
#'   candidate number of subclones.  Defaults to `TRUE`.
#'   
#' @param \dots Further arguments to be passed to [positiveFusedLasso()].
#'   
#' @param verbose A logical indicating whether to print extra information.
#'   Defaults to `FALSE`.
#'   
#' @return A list of k objects of class [c3coFit][c3coFit-class],
#'   where k is the number of candidate number of subclones.
#'   
#' @examples
#' set.seed(7)
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 3
#' bkps <- list(c(100, 250)*10, c(150, 400)*10, c(150, 400)*10)
#' regions <- list(c("(0,3)", "(0,2)", "(1,2)"),
#' c("(1,1)", "(0,1)", "(1,1)"), c("(0,2)", "(0,1)", "(1,1)"))
#' datSubClone <- buildSubclones(len, nbClones, bkps, regions, dataAnnotTP, dataAnnotN)
#' M <- rSparseWeightMatrix(14, 3, sparse.coeff=0.7)
#' dat <- mixSubclones(subClones=datSubClone, M)
#' seg <- segmentData(dat)
#' 
#' l1 <- seq(from=1e-8, to=1e-5, length.out=5)
#' parameters.grid <- list(lambda1=l1,  nb.arch=2:6)
#' fitList <- fitC3co(t(seg$Y1), t(seg$Y2), parameters.grid=parameters.grid)
#' fitListC <- fitC3co(t(seg$Y), parameters.grid=parameters.grid)
#' 
#' 
#' ## A simpler example with toy data
#' l1 <- 1e-4
#' candP <- 2:10
#' parameters.grid <- list(lambda=l1, nb.arch=candP)
#' 
#' dat <- getToyData(n=20, len=100, nbClones=2, nbBkps=5, eps=0.2)  ## almost noiseless!
#' sdat <- dat$segment
#' 
#' res <- fitC3co(sdat$Y, parameters.grid=parameters.grid)
#' pvePlot2(res$config$best, ylim=c(0.6, 1))
#' mus <- sapply(res$fit, FUN=function(x) x@mu)
#' ## the mus are all identical??
#' 
#' @importFrom methods slot
#' @importFrom matrixStats colMaxs
#' @export
fitC3co <- function(Y1, Y2=NULL, parameters.grid=NULL, warn=TRUE, ..., verbose=FALSE) {

    ## preparing data
    Y <- list(Y1=Y1)
    if (!is.null(Y2)) Y$Y2 <- Y2
    ## Sanity checks
    stop_if_not(length(unique(lapply(Y, FUN = dim))) == 1) # are all the dimension equal?

    if (verbose) mprintf("fitC3co() ...\n")
  
    ## problem dimension
    n <- nrow(Y1)
    nseg <- ncol(Y1)

    ## centered version of the data
    Yc <- lapply(Y, FUN = function(y) {
      sweep(y, MARGIN = 1L, STATS = rowMeans(y), FUN = `-`)
    })
    ## Compute total sum of squares
    totVar <- sum((Reduce(`+`, Yc))^2)
    
    ### Define grids
    parameters <- checkParams(parameters.grid, length(Y), nseg, verbose)
    configs <- parameters$configs
    nb.arch <- parameters$nb.arch
    
    fitList <- allRes <- allLoss <- list()
    bestConfigp <- allConfig <- NULL
    for (it in seq_along(nb.arch)) {
        pp <- nb.arch[it]
        if (verbose) mprintf(" - Iteration #%d (%d latent features) of %d ...\n", it, pp, length(nb.arch))

        ## Initialization
        bestRes <- NULL
        bestBIC <- -Inf
        bestConfig <- aConf <- NULL
        ## Z0 are initialized with the centered version of the data
        Z0 <- initializeZ(Yc$Y1, Yc$Y2, p=pp, ...)
        if (verbose) {
            mprintf("   + Parameter configuration: (%s)\n",
                    comma(colnames(configs)))
        }
        allRes[[pp]] <- list()
        allLoss[[pp]] <- list()
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
            res <- positiveFusedLasso(Y = Y, Z = Z0, lambda = cfg)
            
            if (verbose) mprintf(", BIC = ")
            stats <- modelFitStatistics(Reduce(`+`, Y), res@E$Y, res@W, res@S$Z)
            BIC <- stats[["BIC"]]
            if (verbose) mprintf("%g", BIC)
            aConf <- c(pp, cfg, stats[["PVE"]], BIC, stats[["logLik"]], stats[["loss"]])
            ## replace above line by: aConf <- stats
            allConfig <- rbind(allConfig, aConf)
            allRes[[pp]][[cc]] <- res
            allLoss[[pp]][[cc]] <- stats[["loss"]]
            
            if (BIC > bestBIC) { ## BIC has improved: update best model
                bestRes <- res
                bestBIC <- BIC
                bestConfig <- aConf
                if (verbose) mprintf("*")
            }
            if (verbose) mprintf("\n")
        }

        fitList[[it]] <- bestRes
        bestConfigp <- rbind(bestConfigp, bestConfig)
        ## sanity check: minor CN < major CN in the best parameter
        ## configurations (not for all configs by default)
        if (!is.null(Y2) && warn) {
            Z <- slot(bestRes, "S")
            dZ <- Z$Z2 - Z$Z1
            tol <- 1e-2  ## arbitrary tolerance...
            if (min(dZ) < -tol) {
                warning("For model with ", pp, " features, some components in minor latent profiles are larger than matched components in major latent profiles")
            }
        }
	
        if (verbose) mprintf(" - Iteration #%d (%d latent features) of %d ... DONE\n", it, pp, length(nb.arch))
    } ## for (it ...)
    
    names(fitList) <- nb.arch
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


checkParams <- function(parameters.grid, M, nseg, verbose) {
  
    lambda <-  parameters.grid$lambda
    lambda1 <- parameters.grid$lambda1
    lambda2 <- parameters.grid$lambda2
    if (!is.null(lambda)) {
      lambda1 <- lambda
      if (M > 1) {
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
      if (M == 2) {
        if (is.null(lambda1) && is.null(lambda2)) {
          lambda <- 10^(-seq(from=6, to=4, length.out=10))
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
          lambda <- 10^(-seq(from=6, to=4, length.out=10))
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
        nb.arch <- seq(from=2, to=nseg-1, by=1)
        if (verbose) {
            message("Parameter 'nb.arch' not provided. Using default value: ", comma(nb.arch))
        }
    }
  
    list(nb.arch = nb.arch, configs = configs)
}
