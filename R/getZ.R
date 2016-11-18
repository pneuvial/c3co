get.Z <- function(W, Y1, Y2=NULL, lambda1, lambda2) {
  L <- ncol(Y1); p <- ncol(W)
  Dm1 <- -Matrix::bandSparse(L,L,k=0:(L-1))
  X <- kronecker(Dm1,W)

  z1.tilde <- glmnet::glmnet(X,c(Y1), lambda=lambda1, intercept=F, standardize = FALSE)$beta
  Z1 <- apply(matrix(z1.tilde,L,p,byrow=TRUE), 2, function(z) -rev(cumsum(rev(z))))
  if(!is.null(Y2)){
    z2.tilde <- glmnet::glmnet(X,c(Y2), lambda=lambda2, intercept=F, standardize = FALSE)$beta
    Z2 <- apply(matrix(z2.tilde,L,p,byrow=TRUE), 2, function(z) -rev(cumsum(rev(z))))
  }else{
    Z2 <- NULL
  }
  return(list(Z1=Z1,Z2=Z2))
}
