#' @importFrom limSolve lsei
get.W <- function(Zbar, Ybar) {
  p <- ncol(Zbar)
  round(t(apply(Ybar, MARGIN=1L, FUN=function(ybar) {
    lsei(A=Zbar, B=ybar, E=rbind(rep(1, times=p)), F=1,
         H=rep(0, times=p), G=diag(p), tol=1e-15)$X
  })),2)
}
