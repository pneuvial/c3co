get.W <- function(Zbar, Ybar) {
  p <- ncol(Zbar)
  return(t(apply(Ybar, 1, function(ybar) {
    lsei(A = Zbar, B = ybar, E = rbind(rep(1,p)) , F = 1,
         H = rep(0,p), G = diag(p))$X
  })))
}
