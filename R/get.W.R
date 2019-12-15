#' Get the matrix of weights W for fixed values of the subclones Z and data Y
#' 
#' Optimisation is done by solving a least square problem with inequality constraint thanks to \pkg{lsei} package
#' 
#' @param Zt The transpose of Z where Z is a matrix with K rows (number of subclones) and M x J columns (number of signal x number of segments) 
#' 
#' @param Y a matrix with n rows (number of samples) and M x J columns (number of signal x number of segments) 
#'
#' @param type integer code determining algorithm to use 1=lsei, 2=solve.QP from R-package \pkg{quadprog}
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
#' @importFrom limSolve lsei
get.W <- function(Zt, Y, type = 1L) { # partial fix to #58
  K <- ncol(Zt)
  n <- nrow(Y)
  
  ## Sanity checks
  stop_if_not(ncol(Y) == nrow(Zt))
  stop_if_not(K <= ncol(Y))
  
  E <- matrix(rep(1, times = K), nrow = 1L, ncol = K)
  H <- double(K)
  G <- diag(K)

  converged <- TRUE # check if the problem is solvable
  W <- apply(Y, MARGIN = 1L, FUN = function(y) {
    lsei_out <- lsei(A = Zt, B = y, E = E, F = 1, H = H, G = G, type = type)
    converged <<- !lsei_out$IsError
    lsei_out$X
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
