#' @importFrom glmnet glmnet
#' @importFrom Matrix bandSparse
#' @importFrom matrixStats colCumsums
#' @importFrom methods as
get.Z <- function(Y, lambda, W, WtWm1) {
  L <- ncol(Y)
  p <- ncol(W)
  ## temp variables
  X1 <- bandSparse(L, L-1, k=-(1:(L-1)))
  W.WtWm1 <- W %*% WtWm1
  Pw <- W.WtWm1 %*% t(W)
  
  ## Lasso regression 
  X.tilde <- kronecker(scale(X1, center = TRUE, scale = FALSE), as(W, "sparseMatrix"))
  ## Why as.numeric() below? Is it of a different type? /HB 2018-02-27
  y.tilde <- as.numeric(sweep(Y, MARGIN = 1L, STATS = rowMeans(Pw %*% Y), FUN = `-`)) 
  z.tilde <- glmnet(X.tilde, y.tilde, lambda=lambda, intercept=FALSE, standardize=FALSE)$beta
  
  ## Go back to Z 
  Z <- matrix(z.tilde, nrow=L-1, ncol=p, byrow=TRUE)
  X1.Z <- rbind(0, colCumsums(Z))
  
  sweep(X1.Z, MARGIN = 2L, STATS = colMeans(t(Y) %*% W.WtWm1 - X1.Z), FUN = `+`)
}
