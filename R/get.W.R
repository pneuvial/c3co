#' Get the matrix of weights W for fixed values of the subclones Z and data Y
#' 
#' Optimisation is done by solving a least square problem with inequality constraint thanks to \pkg{lsei} package
#' 
#' @param Zt The transpose of Z where Z is a matrix with K rows (number of subclones) and J columns (number of segments) 
#' 
#' @param Y a matrix with n rows (number of samples) and J columns (number of segments) 
#'
#' @return a matrix with n rows (number of samples) and K columns (number of subclones)
#' 
#' @examples 
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
#' E <- matrix(rnorm(n*J, sd = 0), nrow = n, ncol = J)
#'
#' Y <- W %*% Z + E
#'
#' What <- c3co:::get.W(t(Z), Y)
#' 
#' @importFrom limSolve lsei
get.W <- function(Zt, Y) {
  K <- ncol(Zt)

  ## Sanity checks
  stop_if_not(ncol(Y) == nrow(Zt))

  W <- t(apply(Y, MARGIN = 1L, FUN = function(y) {
    lsei(
      A = Zt,
      B = y,
      E = matrix(rep(1, times = K), nrow = 1L),
      F = 1,
      H = rep(0, times = K),
      G = diag(K),
      type = 2L)$X
  }))

  ## Sanity checks
  stop_if_not(nrow(W) == nrow(Y), ncol(W) == K)

  W
}
