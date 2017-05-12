#' @importFrom glmnet glmnet
#' @importFrom Matrix bandSparse
#' @importFrom methods as
get.Z <- function(W, Y, lambda) {
    stopifnot(length(lambda) == 1L)  ## sanity check

    L <- ncol(Y)
    p <- ncol(W)
    Dm1 <- -bandSparse(L, L, k=0:(L-1))
    sW <- as(W, "sparseMatrix")  ## to avoid a NOTE due to kronecker's multiple dispatch
    X <- kronecker(Dm1, sW)
    rm(sW)
    
    z.tilde <- glmnet(X, c(Y), lambda=lambda, intercept=FALSE,
                      standardize=FALSE)$beta

    Z <- matrix(z.tilde, nrow=L, ncol=p, byrow=TRUE)
    Z <- apply(Z, MARGIN=2L, FUN=function(z) {
      -rev(cumsum(rev(z)))
    })

    Z
}

#' @importFrom glmnet glmnet
#' @importFrom Matrix bandSparse
#' @importFrom methods as
#' @examples
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 3
#' bkps <- list(c(100, 250)*10, c(150, 400)*10, c(150, 400)*10)
#' regions <- list(c("(0,3)", "(0,1)", "(1,2)"),
#' c("(1,1)", "(0,1)", "(1,1)"), c("(0,2)", "(0,1)", "(1,1)"))
#' datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN,
#'                               nbClones, bkps, regions)
#' W <- rSparseWeightMatrix(15, 3, 0.4)
#' simu <- mixSubclones(subClones=datSubClone, W)
#' seg <- segmentData(simu)
#' Y1 <- seg$Y1
#' Y2 <- seg$Y2
#' W1 <- cbind(W, 1-rowSums(W))
#' v1 <- get.Z.v2(W1, t(Y1), lambda=1e-5)
#' v2 <- get.Z.v2(W1, t(Y2), lambda=1e-5)
get.Z.v2 <- function(W, Y, lambda) {
  stopifnot(length(lambda) == 1L)  ## sanity check
  
  L <- ncol(Y)
  p <- ncol(W)
  
  ## temp variables
  X1 <- bandSparse(L, L-1, k=-(1:(L-1)))
  W.WtWm1 <- W %*% solve(t(W)%*%W)
  Pw <- W.WtWm1 %*% t(W)
  
  ## Lasso regression 
  X.tilde <- kronecker(scale(X1,TRUE,FALSE), as(W, "sparseMatrix")) 
  y.tilde <- as.numeric(sweep(Y,1,rowMeans(Pw %*% Y),"-")) 
  z.tilde <- glmnet(X.tilde, y.tilde, lambda=lambda, intercept=FALSE, standardize=FALSE)$beta
  
  ## Go back to Z 
  Z <- matrix(z.tilde, nrow=(L-1), ncol=p, byrow=TRUE)
  X1.Z  <- rbind(0,apply(Z, 2, cumsum))
  
  return(sweep(X1.Z, 2, colMeans(t(Y) %*% W.WtWm1 - X1.Z), "+"))
}


