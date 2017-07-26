#' @importFrom limSolve lsei
get.W <- function(Zc, Yc) {
  p <- ncol(Zc)
  t(apply(Yc, MARGIN=1L, FUN=function(y) {
    lsei(A=Zc, B=y, E=rbind(rep(1, times=p)), F=1,
         H=rep(0, times=p), G=diag(p), type=2)$X
  }))
}
