#' @importFrom glmnet glmnet
#' @importFrom Matrix bandSparse
#' @importFrom methods as
get.Z <- function(W, Y, lambda) {
    stopifnot(length(lambda) == 1L)  ## sanity check

    L <- ncol(Y)
    p <- ncol(W)
    Dm1 <- -bandSparse(L, L, k=0:(L-1))
    sW <- as(W, "sparseMatrix")  ## to avoid a NOTE due to kronecker's multiple dispatch
    X <- kronecker(Dm1, sW)
    rm(sW)
    
    z.tilde <- glmnet(X, c(Y), lambda=lambda, intercept=FALSE,
                      standardize=FALSE)$beta

    Z <- matrix(z.tilde, nrow=L, ncol=p, byrow=TRUE)
    Z <- apply(Z, MARGIN=2L, FUN=function(z) {
      -rev(cumsum(rev(z)))
    })

    Z
}
