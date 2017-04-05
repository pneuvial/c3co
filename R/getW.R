#' @import limSolve lsei
get.W <- function(Zbar, Ybar) {
  p <- ncol(Zbar)
  t(apply(Ybar, MARGIN=1L, function(ybar) {
    lsei(A=Zbar, B=ybar, E=rbind(rep(1, times=p)), F=1,
         H=rep(0, times=p), G=diag(p))$X
  }))
}
