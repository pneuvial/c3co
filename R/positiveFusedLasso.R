#' Positive fused lasso function
#'
#' @export
#' @param Y1 A matrix containing the segmented minor copy number (\code{n} patients in row and \code{L} segments in columns)
#' @param Y2 A matrix containing the segmented major copy number (\code{n} patients in row and \code{L} segments in columns)
#' @param nb.arch An integer which is the number of archetypes in the model
#' @param lambda1 A real number, the coefficient for the fused penalty for minor copy numbers 
#' @param lambda2 A real number, the coefficient for the fused penalty for major copy numbers 
#' @param init.random if you want to use random initialization set paramater to TRUE
#' @param eps criterion to stop algorithm (when W do not change sqrt(sum((W-W.old)^2)<eps) 
#' @param max.iter maximum number of iterations of the algorithm
#' @param new.getZ TRUE if you want to parallelize inferrence of Minor and Major copy numbers
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
#' res <- positiveFusedLasso(seg$Y1, seg$Y2, nb.arch=4, lambda1=lambda, lambda2=lambda)
#' showPosFused(res)
#' resC <- positiveFusedLasso(seg$Y1+seg$Y2, Y2=NULL, nb.arch=2, lambda1=lambda, lambda2=lambda)
#'
#' @importFrom parallel mclapply
positiveFusedLasso <- function(Y1, Y2, nb.arch, lambda1, lambda2, init.random=FALSE,
                           eps=1e-2, max.iter=50, verbose=FALSE, new.getZ=FALSE) {
    
    ## problem dimensions
    n <- nrow(Y1) # number of individuals
    L <- ncol(Y1) # number of loci/segments
    ## _______________________________________________________
    ## STEP 0: INITIALIZATION
    if (is.null(Y2)){
        Y <- Y1
    } else {
        Y <- Y1 + Y2
    }
    if (!init.random){
        ## initializing Z by clustering on Y
        hc <- stats::hclust(stats::dist(Y),method="ward.D")
        cluster <- stats::cutree(hc, nb.arch)
        ## averaging the Y over the clusters to initialize the archetypes
        Z1.init <- sapply(split(as.data.frame(Y1), cluster), colMeans)
        if (!is.null(Y2)){
            Z2.init <- sapply(split(as.data.frame(Y2), cluster), colMeans)
        } else {
            Z2.init <- NULL
        }
    } else {
        ii <- sample(1:n,nb.arch, replace=FALSE)
        Z1.init <- t(Y1[ii, ])
        Z2.init <- NULL
        if (!is.null(Y2)){
            Z2.init <- t(Y2[ii, ])
        }
    }
    Z <- list(Z1=Z1.init, Z2=Z2.init)
    ## main loop for alternate optimization
    iter <- 0
    cond <- FALSE
    delta <- Inf
    Ymat <- cbind(Y1, Y2)
    while (!cond) {
        iter <- iter + 1
        if (verbose) message("Iteration: ",iter)
        ## __________________________________________________
        ## STEP 1: optimize over W (fixed Z1, Z2)
        Zmat <- rbind(Z$Z1, Z$Z2)
        W <- get.W(Zmat, Ymat)
        ## __________________________________________________
        ## STEP 2: optimize over Z (fixed W)
        if (!new.getZ){
            Z <- get.Z(W, Y1, Y2, lambda1, lambda2)
        } else {
            Z <- mclapply(list(list(Y=Y1, lambda=lambda1),
                               list(Y=Y2, lambda=lambda2)), 
                          get.Z.new, W=W)
            names(Z) <- c("Z1", "Z2")
        }
        
        ## __________________________________________________
        ## STEP 3: check for convergence of the weights
        if (iter>1) {
            delta <- sqrt(sum((W-W.old)^2))
        }
        cond <- (iter > max.iter | delta < eps)
        
        if (verbose) message("delta:", round(delta, 4))
        W.old <- W
    }
    if (verbose) message("DONE!")
    
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
    kZ <- sum(apply(Z, 2, diff) != 0)
    BIC <-  n*L*log(loss) + kZ*log(n*L)
    PVE <- 1-resVar/predVar
    
    if (!is.null(Y2)) {  ## sanity check: minor CN < major CN
        dZ <- Z2-Z1
        tol <- 1e-2  ## arbitrary tolerance...
        if (min(dZ) < - tol) {
            warning("For model with ", nb.arch, " features, some components in minor latent profiles are larger than matched components in major latent profiles")
        }
    }
    S <- list(Z=Z, Z1=Z1, Z2=Z2)
    E <- list(Y1=Y1hat, Y2=Y2hat)
    param <- list(nb.arch=nb.arch, lambda1=lambda1, lambda2=lambda2)
    objRes <- methods::new("posFused", S=S, W=W, E=E, BIC=BIC, PVE=PVE, param=param)
    return(objRes)
}
