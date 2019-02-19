#' Create a toy subclone and mixture data set
#' 
#' @param n The number of observations.
#'   
#' @param len The number of loci in each subclone.
#'   
#' @param nbClones The number of subclones.
#'   
#' @param nbBkps The total number of breakpoints.
#'   
#' @param eps A numeric value, the signal to noise ratio for simulated data.
#'   
#' @param weightSparsity A numeric value in \eqn{[0,1]}: weights under 
#'   `weightSparsity` are set to 0. This parameter controls the sparsity of
#'   the weight matrix.
#'   
#' @param intercept A logical value indicating whether an intercept should be
#'   added (this corresponds to the presence of a normal subclone).
#' 
#' @return A list with the following elements: 
#'  \describe{
#'  \item{W}{A `n`-by-`nbClones` matrix of weights}
#'  \item{segment}{A list of two elements:
#'   \describe{
#'     \item{Y}{An `n`-by-`nbBkp+1` matrix of observed CN signals}
#'     \item{Z}{An `nbClones`-by-`nbBkp+1` matrix of latent features
#'              (subclones)}}
#'  }
#'  \item{locus}{only returned if \code{returnLocus} is TRUE: A list of two elements:
#'   \describe{
#'     \item{Y}{An `n`-by-`len` matrix of observed CN signals}
#'     \item{Z}{An `nbClones`-by-`len` matrix of latent features
#'              (subclones)}}
#'  }
#'  }
#'   
#' @details For simplicity, the breakpoints positions are drawn uniformly from 
#'   the set of all possible positions.
#'   
#' @export
#' @examples
#' 
#' len <- 100
#' nbClones <- 2
#' nbBkps <- 5
#' n <- 4
#' 
#' dat <- getToyData(n, len, nbClones, nbBkps, eps=0)  ## noiseless
#' matplot(t(dat$locus$Y), t="s")
#' matplot(t(dat$segment$Y), t="s")
#' 
#' dat <- getToyData(n, len, nbClones, nbBkps, eps=0.2)  ## noisy
#' matplot(t(dat$locus$Y), t="l")
#' matplot(t(dat$segment$Y), t="s")
#' 
#' len <- 1000
#' nbClones <- 5
#' nbBkps <- 10
#' eps <- 1
#' n <- 20
#' 
#' dat <- getToyData(n, len, nbClones, nbBkps, eps=0)  ## noiseless
#' matplot(t(dat$locus$Y), t="s")
#' matplot(t(dat$segment$Y), t="s")
#' 
#' dat <- getToyData(n, len, nbClones, nbBkps, eps=0.2)  ## noisy
#' matplot(t(dat$locus$Y), t="l")
#' matplot(t(dat$segment$Y), t="s")
#' 
#' l1 <- seq(from=1e-6, to=1e-4, length.out=10)
#' parameters.grid <- list(lambda=l1, nb.arch=2:9)
#' 
#' Y <- dat$segment$Y
#' fit <- fitC3co(Y, parameters.grid=parameters.grid)
#' pvePlot2(fit$config$best)
#' 
getToyData <- function(n, len, nbClones, nbBkps, eps, weightSparsity = 0.1, intercept = TRUE, returnLocus = TRUE) {
    
    ## number of segments
    nbSegs <- nbBkps + 1 
    
    ## breakpoint positions
    bkp <- sample(len, size = nbBkps, replace = FALSE)
    bkp <- sort(bkp)
    
    ## segment lengths
    segLens <- diff(c(0, bkp))
    segLens <- c(segLens, len - sum(segLens))
    
    if (intercept) {    
        nbClones <- nbClones+1
    }
    ## weights
    ru <- runif(n*nbClones)
    ru[ru < weightSparsity] <- 0
    W <- matrix(ru, nrow = n, ncol = nbClones)
    W <- round(sweep(W, MARGIN = 1L, STATS = rowSums(W), FUN = `/`), digits = 2L)
    W[, nbClones] <- 1 - rowSums(W[, -nbClones, drop = FALSE]) ## make sure rows sum to 1 after rounding
    
    Zlist <- list()
    ## subclones (segment-level)
    z <- rnorm(nbSegs*nbClones)
    Zs <- round(matrix(z, nrow = nbClones, ncol = nbSegs))
    if (intercept) {  ## last clone should be 'normal'
        Zs[nbClones, ] <- 1
    }
    
    ## subclones (locus-level)
    Zl <- t(apply(Zs, MARGIN = 1L, FUN = rep, times = segLens))

    ## segment-level data
    e <- rnorm(n*nbSegs, sd=eps)
    Es <- matrix(e, nrow=n, ncol=nbSegs)  ## noise
    Ys <- W %*% Zs + Es                   ## observations
    seg <- list(Y=Ys, Z=Zs)
    
    res <- list(W = W, segment = seg)

    ## locus-level data
    if (returnLocus) {
        e <- rnorm(n*len, sd = eps)
        El <- matrix(e, nrow = n, ncol = len)  ## noise
        Yl <- W %*% Zl + El                ## observations
        loc <- list(Y = Yl, Z = Zl)
        res$locus <- loc
    }
    
    res
}
