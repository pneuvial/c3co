#' Positive fused lasso function
#'
#' @export
#' @param Y1 A matrix containing the segmented minor copy number (\code{n} patients in row and \code{L} segments in columns)
#' @param Y2 A matrix containing the segmented major copy number (\code{n} patients in row and \code{L} segments in columns)
#' @param nb.arch An integer which is the number of archetypes in the model
#' @param lambda1 A real number which is the penalty coefficient for the fused lasso on the minor copy number dimension
#' @param lambda2 A real number which is the penalty coefficient for the fused lasso on the major copy number dimension
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
#' M = matrix(c(40,30, 0,0,70, 15,10, 0, 35,15,0,0 ,0,0,0), byrow=TRUE, ncol=3)
#' simu <- mixSubclones(subClones=datSubClone,M)
#' YTCNtoSeg <- t(sapply(simu, function(cc) cc$tcn))
#' YDHtoSeg <- t(sapply(simu, function(cc) cc$dh))
#' sim <- t(rbind(YTCNtoSeg,YDHtoSeg))
#' res <- jointseg::jointSeg(sim, method = "RBS", K=8)
#' bkp <- res$bestBkp
#' Y1 <- t(sapply(simu, function(cc) cc$c1))
#' Y2 <-  t(sapply(simu, function(cc) cc$c2))
#' Y1seg <- t(apply(Y1, 1, matrixStats::binMeans, x=1:ncol(Y1),bx= c(1,bkp,ncol(Y1)), na.rm=TRUE))
#' Y2seg <- t(apply(Y2, 1, matrixStats::binMeans, x=1:ncol(Y1),bx= c(1,bkp,ncol(Y1)), na.rm=TRUE))
#' lambda <- 0.001
#' system.time(rC1C2 <- positive.fused(Y1seg,Y2seg, 4,lambda1 = lambda, lambda2 = lambda))
#' system.time(rC1C2new <- positive.fused(Y1seg,Y2seg, 4,
#' lambda1 = lambda, lambda2 = lambda, new.getZ=TRUE))
#' rTCN <- positive.fused(Y1seg+Y2seg,NULL, 2,lambda1 = lambda, lambda2 = lambda)
#' showPosFused(rTCN)
positive.fused <- function(Y1, Y2, nb.arch, lambda1, lambda2, init.random=FALSE,
                           eps = 1e-2, max.iter = 50, verbose=F, new.getZ=FALSE) {
  
  ## problem dimensions
  n <- nrow(Y1) # number of individuals
  L <- ncol(Y1) # number of loci
  ## _______________________________________________________
  ## STEP 0: INITIALIZATION
  ## Under recommandations of NMF and Archetypal analysis random initialization to find the global minimum.
  if(is.null(Y2)){
    Y=Y1
  }else{
    Y=Y1+Y2
  }
  if(!init.random){
    ## initializing Z by clustering on Y
    cluster <- stats::cutree(stats::hclust(stats::dist(Y),method="ward.D"), nb.arch)
    ## averaging the Y over the clusters to initialize the archetypes
    Z1.init <- sapply(split(as.data.frame(Y1), cluster), colMeans)
    if(!is.null(Y2)){
      Z2.init <- sapply(split(as.data.frame(Y2), cluster), colMeans)
    }else{
      Z2.init <- NULL
    }
  }else{
    i <- sample(1:n,nb.arch, replace=FALSE)
    Z1.init <- t(Y1[i,])
    if(!is.null(Y2)){
      Z2.init <- t(Y2[i,])
    }else{
      Z2.init <- NULL
    }
  }
  Z <- list(Z1 = Z1.init, Z2 = Z2.init)
  ## main loop for alternate optimization
  iter <- 0
  cond <- FALSE
  delta <- Inf
  min.loss <- 1e5
  while(!cond) {
    iter <- iter + 1
    if(verbose) cat("\niter number ",iter)
    ## __________________________________________________
    ## STEP 1: optimize over W (fixed Z1, Z2)
    W <- get.W(rbind(Z$Z1,Z$Z2), cbind(Y1,Y2))
    ## __________________________________________________
    ## STEP 2: optimize over Z (fixed W)
    if(!new.getZ){
      Z <- get.Z(W, Y1, Y2, lambda1, lambda2)
    }else{
      ## test
      Z <- parallel::mclapply(list(list(Y=Y1,lambda=lambda1),list(Y=Y2,lambda=lambda2)), get.Z.new, W=W)
      names(Z)= c("Z1", "Z2")
    }
    
    ## __________________________________________________
    ## STEP 3: check for convergence of the weights
    if (iter>1) {
      delta <- sqrt(sum((W-W.old)^2))}
    
    cond <- (iter > max.iter | delta < eps)
    
    if(verbose) cat("\ndelta =",round(delta, 4))
    W.old <- W
  }
  if(!is.null(Y2)){
    loss <- (sum((Y1-W %*% t(Z$Z1))^2)+sum((Y2-W %*% t(Z$Z2))^2))/(n*L)
    loss <- (sum((Y1+Y2-W %*% t(Z$Z1+Z$Z2))^2))/(n*L)
    PVE <- 1-(sum((Y1+Y2-W %*% t(Z$Z1+Z$Z2))^2))/(sum(((Y1+Y2)-rowMeans(Y1+Y2))^2))
  }else{
    loss <- (sum((Y1-W %*% t(Z$Z1))^2))/(n*L)
    PVE <- 1-(sum((Y1-W %*% t(Z$Z1))^2))/(sum((Y1-rowMeans(Y1))^2))
  }
  if(loss<min.loss){
    if(!is.null(Y2)){
      res <- list(Z=Z$Z1+Z$Z2, Z1=Z$Z1, Z2=Z$Z2, W=W, Y.hat=list(Y1 = W %*% t(Z$Z1), Y2 = W %*% t(Z$Z2)))
    }else{
      res <- list(Z=Z$Z1, Z1=Z$Z1, Z2=Z$Z2, W=W, Y.hat=list(Y1 = W %*% t(Z$Z1), Y2 = NULL))
    }
    min.loss <- loss
  }
  if(verbose) cat("\nDONE!\n")
  kZ <- sum(apply(res$Z, 2, diff)!=0)
  BIC <-  n*L*log(loss)+kZ*log(n*L)
  
  idx <- ncol(res$Z2)
  checkC2supC1 <- sapply(idx, function(ii){
    jj <- which(res$Z1[,ii]>res$Z2[,ii])
    if(length(jj)>0){
      warning(sprintf("For model with %s features, some components in minor latent profiles are larger than matched components in major latent profiles", nb.arch))
    }
  })
  
  objRes <- methods::new("posFused", S = list(Z = res$Z, Z1 = res$Z1, Z2 = res$Z2), W = res$W, E = list(Y1 = res$Y.hat$Y1, Y2 = res$Y.hat$Y2), BIC=BIC, PVE=PVE, param=list(nb.arch=nb.arch, lambda1=lambda1, lambda2=lambda2))
  return(objRes)
  
}


