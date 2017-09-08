#' @importFrom glmnet glmnet
#' @importFrom Matrix bandSparse
#' @importFrom methods as
get.Z <- function(Y, W, WtWm1, lambda) {
  
  L <- ncol(Y)
  p <- ncol(W)
  ## temp variables
  X1 <- bandSparse(L, L-1, k=-(1:(L-1)))
  W.WtWm1 <- W %*% WtWm1
  Pw <- W.WtWm1 %*% t(W)
  
  ## Lasso regression 
  X.tilde <- kronecker(scale(X1,TRUE,FALSE), as(W, "sparseMatrix")) 
  y.tilde <- as.numeric(sweep(Y, 1, rowMeans(Pw %*% Y), "-")) 
  z.tilde <- glmnet(X.tilde, y.tilde, lambda=lambda, intercept=FALSE, standardize=FALSE)$beta
  
  ## Go back to Z 
  Z <- matrix(z.tilde, nrow=(L-1), ncol=p, byrow=TRUE)
  X1.Z  <- rbind(0, apply(Z, 2, cumsum))
  
  return(sweep(X1.Z, 2, colMeans(t(Y) %*% W.WtWm1 - X1.Z), "+"))
}
