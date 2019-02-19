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
#' nbSegments <- 11
#' nbClones   <- 4
#' 
#' Z <- matrix(1, nrow = nbClones, ncol = nbSegments)
#' Z[2, 2] <- 2
#' Z[3, 5:6] <- 2
#' Z[4, 9:10] <- 2
#'
#' W <- diag(rep(1, times = nrow(Z)))
#' E <- matrix(rnorm(nrow(W)*ncol(Z), sd = 0), nrow = nrow(W), ncol = ncol(Z))
#' Y <- W %*% Z + E
#' W <- c3co:::get.W(t(Z), Y)
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
