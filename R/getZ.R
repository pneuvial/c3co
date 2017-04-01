get.Z <- function(W, Y, lambda) {
    stopifnot(length(lambda)==1L)  ## sanity check
    
    L <- ncol(Y);
    p <- ncol(W)
    Dm1 <- -Matrix::bandSparse(L, L, k=0:(L-1))
    X <- kronecker(Dm1, W)

    z.tilde <- glmnet::glmnet(X, c(Y), lambda=lambda, intercept=FALSE, standardize=FALSE)$beta
    Z <- apply(matrix(z.tilde, L, p, byrow=TRUE), 2, function(z) -rev(cumsum(rev(z))))
    return(Z)
}
