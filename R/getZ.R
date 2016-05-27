get.Z <- function(W, Y1, Y2=NULL, lambda1, lambda2) {
  L <- ncol(Y1); p <- ncol(W)
  Dm1 <- -bandSparse(L,L,k=0:(L-1))
  X <- kronecker(Dm1,W)

  z1.tilde <- glmnet(X,c(Y1), lambda=lambda1, intercept=F, standardize = FALSE)$beta
  Z1 <- apply(matrix(z1.tilde,L,p,byrow=TRUE), 2, function(z) -rev(cumsum(rev(z))))
  ## sanity check : check if z1>0
  Z1 <- apply(Z1, 2, function(z){
   ii <- which(z<0)
   z[ii] <- 0
   return(z)
  })
  if(!is.null(Y2)){
    z2.tilde <- glmnet(X,c(Y2), lambda=lambda2, intercept=F, standardize = FALSE)$beta
    Z2 <- apply(matrix(z2.tilde,L,p,byrow=TRUE), 2, function(z) -rev(cumsum(rev(z))))
    ## sanity check : check if (z2>z1)
    idx <- 1 : ncol(Z2)
    Z1 <- sapply(idx, function(ii){
      jj <- which(Z1[,ii]>Z2[,ii])
      Z1[jj, ii] <- Z2[jj, ii]
      return(Z1[,ii])
    })
  }else{
    Z2 <- NULL
  }

  return(list(Z1=Z1,Z2=Z2))
}
