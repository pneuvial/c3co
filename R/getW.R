#' @importFrom limSolve lsei
get.W <- function(Zc, Yc) {
  p <- ncol(Zc)
  E <- matrix(rep(1, times = p), nrow = 1L)
  H <- rep(0, times = p)
  G <- diag(p)
  t(apply(Yc, MARGIN = 1L, FUN=function(y) {
    lsei(A = Zc, B = y, E = E, F = 1, H = H, G = G, type = 2L)$X
  }))
}
