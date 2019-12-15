#' Get the matrix of weights W for fixed values of the subclones Z and data Y
#' 
#' Optimisation is done by solving a least square problem with equality + inequality constraint thanks to the \pkg{quadprog} package
#' 
#' @param Zt The transpose of Z where Z is a matrix with K rows (number of subclones) and M x J columns (number of signal x number of segments) 
#' 
#' @param Y a matrix with n rows (number of samples) and M x J columns (number of signal x number of segments) 
#' 
#' @param nu positive double standing for ridge regularization in the quadratic programming solver. Default is 1e-8.
#'
#' @return a matrix with n rows (number of samples) and K columns (number of subclones)
#' 
#' @examples 
#' ## Example with M = 1 signal
#' J <- 11L  ## Number of segments
#' K <- 4L   ## Number of subclones
#' n <- 4L   ## Number of samples
#' 
#' Z <- matrix(1, nrow = K, ncol = J)
#' Z[2,    2] <- 2
#' Z[3,  5:6] <- 2
#' Z[4, 9:10] <- 2
#'
#' W <- matrix(0, nrow = n, ncol = K)
#' W[1,1] <- W[2,2] <- W[3,3] <- W[4,4] <- 1
#'
#' E <- matrix(rnorm(n*J, sd = 0.1), nrow = n, ncol = J)
#'
#' Y <- W %*% Z + E
#'
#' What <- c3co:::get.W(t(Z), Y)
#' 
#' @importFrom quadprog solve.QP
get.W <- function(Zt, Y, nu = 1e-8) {
  K <- ncol(Zt)
  n <- nrow(Y)
  
  ## Sanity checks
  stop_if_not(ncol(Y) == nrow(Zt))
  stop_if_not(K <= ncol(Y))

  ## Definition od the quadratic program
  Dmat <- crossprod(Zt) + diag(nu, ncol = K, nrow = K)
  Amat <- cbind(rep(1L, times = K), diag(1L, ncol = K, nrow = K))
  bvec <- c(1.0, double(K))

  W <- apply(Y, MARGIN = 1L, FUN = function(y) {
      dvec <- crossprod(Zt, y)
      w <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)$solution
      w[abs(w) < sqrt(.Machine$double.eps)] <- 0
      w
  })

  if (!converged) return(matrix(NA_real_, nrow=n, ncol=K))

  ## WORKAROUND: https://github.com/pneuvial/c3co/issues/52
  if (K == 1L) {
    dim(W) <- c(n, K)
  } else {
    W <- t(W)
  }

  ## Sanity checks
  stop_if_not(nrow(W) == n, ncol(W) == K)

  W
}
