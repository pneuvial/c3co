#' @importFrom glmnet glmnet
#' @importFrom Matrix bandSparse
get.Z <- function(W, Y, lambda) {
    stopifnot(length(lambda) == 1L)  ## sanity check
    
    L <- ncol(Y)
    p <- ncol(W)
    Dm1 <- -bandSparse(L, L, k=0:(L-1))
    X <- kronecker(Dm1, W)

    z.tilde <- glmnet(X, c(Y), lambda=lambda, intercept=FALSE,
                      standardize=FALSE)$beta

    Z <- matrix(z.tilde, nrow=L, ncol=p, byrow=TRUE)
    Z <- apply(Z, MARGIN=2L, FUN=function(z) {
      -rev(cumsum(rev(z)))
    })
    
    Z
}
