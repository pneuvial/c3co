#' Generate weight matrix
#' 
#' @param nb.samp An integer that is the number of samples (number of rows).
#'
#' @param nb.arch An integer that is the number of archetypes
#' (number of columns).
#'
#' @param sparse.coeff A numeric in \eqn{[0,1]} that control the sparsity by
#' rows (the proportion of non-zero entries among all matrix entries).
#'
#' @return A matrix of weights (each row sums to one).
#' 
#' @examples
#' W <- rSparseWeightMatrix(nb.samp=10L, nb.arch=7L, sparse.coeff=0.5)
#'
#' @importFrom Matrix rsparsematrix 
#' @importFrom Matrix rowSums
#' @importFrom Matrix Matrix 
#' @importFrom stats runif
#' @export
rSparseWeightMatrix <- function(nb.samp, nb.arch, sparse.coeff = max(nb.samp, nb.arch)/(nb.samp*nb.arch)) {
  ## Create a sparse matrix for tumor proportions
  m.tum <- rSpMatrix(nrow = nb.samp, ncol = nb.arch,
                     sparse.coeff = sparse.coeff)
  ## Create the vector contamination by normal cell with a uniform distribution between 0.01 and expand contam.max parameter
  x.contam <- runif(nb.samp, min = 0.01, max = 1)
  m <- cbind(m.tum, x.contam)  ## A sparse matrix
  ## Compute the sum of rows (tum+contam)
  tot <- Matrix::rowSums(m)
  ## Round and divide the matrix by the total remove contamination column
  m.tum.res <- round((m/tot)[,-(nb.arch + 1)], digits = 2L)
  Matrix(m.tum.res, sparse = TRUE)
}

#' Generate sparse matrix with at least one non-zero element by rows and by columns
#' @param nrow An integer that is the number of rows.
#'
#' @param ncol An integer that is the number of columns.
#'
#' @param sparse.coeff A numeric in \eqn{[0,1]} that control the sparsity by
#' rows (the proportion of non-zero entries among all matrix entries).
#' Non-zero entries is define by `nnz = sparse.coeff * nrow * ncol` and
#' need to be larger than `nrow` and `ncol`.
#'
#' @return A matrix of weights (each row sums to one).
#' 
#' @examples
#' M <- rSpMatrix(nrow=10L, ncol=7L, sparse.coeff=0.5)
#'
#' @importFrom stats runif
#' @importFrom Matrix spMatrix 
#' @export
rSpMatrix <- function(nrow, ncol, sparse.coeff) {
  rand.x <- function(n) runif(n, min = 0, max = 1)
  nnz <- ceiling(sparse.coeff*nrow*ncol)
  nnz <- as.integer(nnz)
  stop_if_not(nnz >= 0, nrow >= 0, ncol >= 0, nnz >= max(nrow, ncol))
  spMatrix(nrow, ncol,
           i = c(sample(nrow, size = nrow      , replace = FALSE ),
                 sample(nrow, size = nnz - nrow, replace = TRUE )),
           j = c(sample(ncol, size = nnz - ncol, replace = TRUE  ),
                 sample(ncol, size = ncol      , replace = FALSE)),
           x = rand.x(nnz))
}
