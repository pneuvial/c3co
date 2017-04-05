#' Generate weight matrix
#'
#' @param prop.max An integer that is the maximim proportion for a present
#' archetype
#'
#' @param prob.min An integer that is the minimum proportion  for a present
#' archetype
#'
#' @param nb.arch  An interger that is the number of archetypes
#' (number of columns)
#'
#' @param nb.samp An interger that is the number of samples (number of rows)
#'
#' @param contam.max An interger between 0 and 100 that controls the maximal
#' level of contamination by normal cells
#'
#' @param contam.coeff A numeric between 0 and 1 that control the contamination
#' level
#'
#' @param sparse.coeff A numeric between 0 and 1 that control the sparsity by
#' rows
#'
#' @return A matrix of weights (the sum of rows is equal to 1)
#'
#' @examples
#' M <- getWeightMatrix(70, 20, 7, 30)
#'
#' @importFrom stats rbinom
#' @export
getWeightMatrix <- function(prop.max, prob.min, nb.arch, nb.samp,
                            sparse.coeff=0.5, contam.coeff=0.6, contam.max=30) {
    M <- matrix(0, ncol=nb.arch, nrow=nb.samp)
    nbSparse <- replicate(nb.samp, rbinom(1, size=nb.arch-1, prob=sparse.coeff)+1)
    contamVector <- ceiling(rbinom(nb.samp, size=contam.max, prob=contam.coeff)/5)*5
    for (ss in seq_len(nb.samp)) {
        ind <- sample(size=nbSparse[ss], x=1:nb.arch)
        for (ii in ind) {
            if (sum(M[ss, ]) >= 100) {
                M[ss, ii] <- 100-sum(M[ss, ]-contamVector[ss])
            }
            M[ss, ii] <- ceiling(rbinom(1, size=(100-contamVector[ss])-sum(M[ss, ]), prob=0.65)/5)*5
        }
    }
    return(M)
}
