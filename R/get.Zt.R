#' Get the matrix of subclones Z for fixed values of the weights W
#' 
#' Optimisation is done by solving a L1-penalized least-square problem thanks to the \pkg{glmnet} package
#' 
#' @param Y a matrix with n rows (number of samples) and J columns (number of segments) 
#'
#' @param lambda a positive scalar tuning the penalty level for fusion
#' 
#' @param W a matrix with n rows (number of samples) and K columns (number of subclones)
#' 
#' @param WtWm1 a K x K square matrix (with K the number of subclones), precomputed to save time
#'
#' @return The transpose of Z, where Z a matrix with K rows (number of subclones) and
#' J columns (number of segments) 
#' 
#' @examples
#' J <- 11L  ## Number of segments
#' K <- 4L   ## Number of subclones
#' n <- 4L   ## Number of samples
#' 
#' W <- diag(rep(1, times = nrow(Z)))
#' E <- matrix(rnorm(nrow(W)*ncol(Z), sd = 0.1), nrow = nrow(W), ncol = ncol(Z))
#' WtWm1 <- diag(rep(1, times = nrow(Z)))
#' Z <- matrix(1, nrow = K, ncol = J)
#' Z[2,    2] <- 2
#' Z[3,  5:6] <- 2
#' Z[4, 9:10] <- 2
#'
#' W <- matrix(0, nrow = n, ncol = K)
#' W[1,1] <- W[2,2] <- W[3,3] <- W[4,4] <- 1
#' WtWm1 <- tcrossprod(backsolve(qr.R(qr(W)), x=diag(K)))
#'
#' E <- matrix(rnorm(n*J, sd = 0), nrow = n, ncol = J)
#'
#' Y <- W %*% Z + E
#' 
#' Zthat <- c3co:::get.Zt(Y, lambda = 0.01, W = W, WtWm1 = WtWm1)
#' 
#' @importFrom glmnet glmnet
#' @importFrom Matrix bandSparse
#' @importFrom matrixStats colCumsums
#' @importFrom methods as
get.Zt <- function(Y, lambda, W, WtWm1) {
  L <- ncol(Y)  ## FIXME: Renamed 'L' to 'J'
  p <- ncol(W)  ## FIXME: Renamed 'p' to 'K'

  ## Sanity checks
  stop_if_not(length(lambda) == 1L)  
  stop_if_not(nrow(Y) == nrow(W))
  stop_if_not(nrow(WtWm1) == p, ncol(WtWm1) == p)

  ## temp variables
  X1 <- bandSparse(L, L-1L, k=-(1:(L-1)))
  W.WtWm1 <- W %*% WtWm1
  Pw <- W.WtWm1 %*% t(W)
  
  ## FIXME: include the penalty factor in get.Zt
  ## FIXME to be consistent with what is written,
  ## FIXME lambda should be lambda / (2*n*J) when calling glmnet
  
  ## Lasso regression 
  X.tilde <- kronecker(scale(X1, center = TRUE, scale = FALSE), as(W, "sparseMatrix"))
  ## Why as.numeric() below? Is it of a different type? /HB 2018-02-27
  y.tilde <- as.numeric(sweep(Y, MARGIN = 1L, STATS = rowMeans(Pw %*% Y), FUN = `-`)) 
  z.tilde <- glmnet(X.tilde, y.tilde, lambda = lambda, intercept = FALSE, standardize = FALSE)$beta
  
  ## Go back to Z 
  Z <- matrix(z.tilde, nrow=L-1L, ncol=p, byrow=TRUE)
  X1.Z <- rbind(0, colCumsums(Z))
  
  Zt <- sweep(X1.Z, MARGIN = 2L, STATS = colMeans(t(Y) %*% W.WtWm1 - X1.Z), FUN = `+`)

  ## Sanity checks
  stop_if_not(nrow(Zt) == L, ncol(Zt) == p)

  Zt
}
