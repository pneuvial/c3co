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
#' @example 
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 3
#' bkps <- list(c(100, 250)*10, c(150, 400)*10, c(150, 400)*10)
#' regions <- list(c("(0,3)", "(0,1)", "(1,2)"),
#' c("(1,1)", "(0,1)", "(1,1)"), c("(0,2)", "(0,1)", "(1,1)"))
#' datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN,
#'                               nbClones, bkps, regions)
#' W <- rSparseWeightMatrix(15, 7, 0.4)
#' simu <- mixSubclones(subClones=datSubClone, M)
#' seg <- segmentData(simu)
#' Y1 <- seg$Y1
#' v2 <- get.Z.v2(W, t(Y1), lambda=1e-5)
#' v1 <- c3co:::get.Z(W, t(Y1), lambda=1e-5)



get.Z.v2 <- function(W, Y, lambda) {
  stopifnot(length(lambda) == 1L)  ## sanity check
  
  L <- ncol(Y)
  p <- ncol(W)
  
  ## temp variables
  X1 <-  bandSparse(L, L-1, k=-(1:(L-1)))
  Il <- bandSparse(L, L, k=0)
  Pw <- W%*%t(solve(t(W)%*%W))%*%t(W)
  Pl <- matrix(1, L, L)/L
  sW <- as(W, "sparseMatrix")  ## to avoid a NOTE due to kronecker's multiple dispatch
  
  ## Lasso regression 
  X.tilde <- kronecker((Il-Pl)%*%X1, sW)
  rm(sW)
  Y.tilde <- Y- kronecker(Pl, Pw)%*%Y
  z.tilde <- glmnet(X.tilde, as.numeric(Y.tilde), lambda=lambda, intercept=FALSE,
                    standardize=FALSE)$beta
  
  ## Go back to Z 
  Z <- matrix(z.tilde, nrow=(L-1), ncol=p, byrow=TRUE)
  gamma <- Pl%*%t(Y-W%*%t(X1%*%Z))%*%(W%*%solve(t(W)%*%W))
  Z = gamma+X1%*%Z
  
  Z
}

