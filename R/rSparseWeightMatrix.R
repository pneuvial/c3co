#' Generate weight matrix
#' @param nb.samp An interger that is the number of samples (number of rows)
#'
#' @param nb.arch  An interger that is the number of archetypes
#' (number of columns)
#'
#' @param sparse.coeff A numeric in  [0, 1] that control the sparsity by
#' rows (the proportion of non-zero entries among all matrix entries)
#'
#' @return A matrix of weights (the sum of rows is equal to 1)
#' 
#' @examples
#' M <- rSparseWeightMatrix(nb.samp=10, nb.arch=7, sparse.coeff=0.5)
#'
#' @importFrom Matrix rsparsematrix 
#' @importFrom Matrix Matrix 
#' @importFrom stats runif
#' @export
rSparseWeightMatrix <- function(nb.samp, nb.arch, sparse.coeff= max(nb.samp, nb.arch)/(nb.samp*nb.arch)){
  ## Create a sparse matrix for tumor proportions
  m.tum <- rSpMatrix(nb.samp, nb.arch, sparse.coeff)
  ## Create the vector contamination by normal cell with a uniform distribution between 0.01 and expand contam.max parameter
  x.contam <- runif(nb.samp, 0.01, 1)
  m <- cbind(m.tum, x.contam)
  ## Compute the sum of rows (tum+contam)
  tot <- rowSums(as.matrix(m))
  ## Round and divide the matrix by the total remove contamination column
  m.tum.res <- round((m/tot)[,-(nb.arch+1)],2)
  return(Matrix(m.tum.res, sparse=TRUE))
}

#' Generate sparse matrix with at least one non-zero element by rows and by columns
#' @param nrow An interger that is the number of rows
#'
#' @param ncol  An interger that is the number of columns
#'
#' @param sparse.coeff A numeric in  [0, 1] that control the sparsity by
#' rows (the proportion of non-zero entries among all matrix entries). Non zero entries is define by  \code{nnz}=\code{sparse.coeff*nrow*ncol} and need to be larger than \code{nrow} and \code{ncol}
#'
#' @return A matrix of weights (the sum of rows is equal to 1)
#' 
#' @examples
#' M <- rSpMatrix(10, 7, 0.5)
#'
#' @importFrom stats runif
#' @importFrom Matrix spMatrix 
#' @export
rSpMatrix <- function(nrow, ncol, sparse.coeff)
{
  nnz <- ceiling(sparse.coeff*nrow*ncol)
  rand.x = function(n) runif(n, 0, 1)
  stopifnot((nnz <- as.integer(nnz)) >= 0,
            nrow >= 0, ncol >= 0,nnz >= max(nrow, ncol))
  spMatrix(nrow, ncol,
           i = c(sample(nrow, nrow, replace = FALSE),sample(nrow, nnz-nrow, replace=TRUE)),
           j = c(sample(ncol, (nnz-ncol), replace = TRUE),sample(ncol, ncol, replace = FALSE)),
           x = rand.x(nnz))
}