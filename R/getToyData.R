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
#' @param dimension An integer value in {1,2}, the dimension of the signals to
#'   be generated (eg 1 for total copy numbers and 2 for minor and major copy
#'   numbers)
#' 
#'   
#' @param intercept A logical value indicating whether an intercept should be
#'   added (this corresponds to the presence of a normal subclone).
#'   
#' @param returnLocus A logical value indicating whether the locus-level data
#'   should be returned. Defaults to TRUE
#'   
#' @return 
#'  \item{W}{A `n`-by-`nbClones` matrix of weights}
#'  \item{segment}{A list of two elements:
#'   \describe{
#'     \item{Y}{An `n`-by-`nbBkp+1` matrix of observed CN signals if \code{dimension==1}, or a list of two such matrices if \code{dimension==2}}
#'     \item{Z}{An `nbClones`-by-`nbBkp+1` matrix of latent features
#'              (subclones)} if \code{dimension==1}, or a list of two such matrices if \code{dimension==2}}
#'  }
#'  \item{locus}{only returned if \code{returnLocus} is TRUE: A list of two elements:
#'   \describe{
#'     \item{Y}{An `n`-by-`len` matrix of observed CN signals, if \code{dimension==1}, or a list of two such matrices if \code{dimension==2}}
#'     \item{Z}{An `nbClones`-by-`len` matrix of latent features
#'              (subclones)}, if \code{dimension==1}, or a list of two such matrices if \code{dimension==2}}
#'  }
#'   
#' @details For simplicity, the breakpoints positions are drawn uniformly from 
#'   the set of all possible positions.
#'   
#' @export
#' @examples
#' 
#' len <- 100L
#' nbClones <- 3L
#' nbBkps <- 5L
#' n <- 10L
#' 
#' dat <- getToyData(n, len, nbClones, nbBkps, eps = 0)  ## noiseless
#' matplot(t(dat$locus$Y), t = "s")
#' matplot(t(dat$segment$Y), t = "s")
#' 
#' dat <- getToyData(n, len, nbClones, nbBkps, eps = 0.2)  ## noisy
#' matplot(t(dat$locus$Y), t="s")
#' matplot(t(dat$segment$Y), t="s")
#' 
#' 
#' \dontrun{
#' l1 <- seq(from = 1e-6, to = 1e-4, length.out = 10L)
#' parameters.grid <- list(lambda = l1, nb.arch = 2:6)
#' Y <- dat$segment$Y
#' fit <- fitC3co(Y, parameters.grid=parameters.grid)
#' pvePlot2(fit$config$best)
#' }
getToyData <- function(n, len, nbClones, nbBkps, eps, 
                       weightSparsity = 0.1, dimension = 1L,
                       intercept = TRUE, returnLocus = TRUE) {
    ## sanity check
    stop_if_not(dimension %in% 1:2)
    stop_if_not(weightSparsity >= 0, weightSparsity <= 1)
    
    ## number of segments
    nbSegs <- nbBkps + 1L
    
    ## breakpoint positions
    bkp <- sample(len, size = nbBkps, replace = FALSE)
    bkp <- sort(bkp)
    
    ## segment lengths
    segLens <- diff(c(0, bkp))
    segLens <- c(segLens, len - sum(segLens))
    
    if (intercept) {    
        nbClones <- nbClones + 1L
    }
    ## weights
    ru <- runif(n*nbClones)
    ru[ru < weightSparsity] <- 0
    W <- matrix(ru, nrow = n, ncol = nbClones)
    W <- round(sweep(W, MARGIN = 1L, STATS = rowSums(W), FUN = `/`), digits = 2L)
    W[, nbClones] <- 1 - rowSums(W[, -nbClones, drop = FALSE]) ## make sure rows sum to 1 after rounding
    
    YsegList <- list()
    YlocList <- list()
    ZsegList <- list()
    ZlocList <- list()
    for (dd in 1:dimension) {
        ## subclones (segment-level)
        z <- round(rnorm(nbSegs*nbClones))
        Zs <- matrix(z, nrow = nbClones, ncol = nbSegs)
        if (intercept) {  ## last clone should be 'normal'
            Zs[nbClones, ] <- 1
        }
        
        ## subclones (locus-level)
        Zl <- t(apply(Zs, MARGIN = 1L, FUN = rep, times = segLens))

        ## segment-level data
        e <- rnorm(n*nbSegs, sd=eps)
        Es <- matrix(e, nrow=n, ncol=nbSegs)  ## noise
        Ys <- W %*% Zs + Es                   ## observations
        
        YsegList[[dd]] <- Ys
        ZsegList[[dd]] <- Zs

        ## locus-level data
        if (returnLocus) {
            e <- rnorm(n*len, sd = eps)
            El <- matrix(e, nrow = n, ncol = len)  ## noise
            Yl <- W %*% Zl + El                ## observations
            
            YlocList[[dd]] <- Yl
            ZlocList[[dd]] <- Zl
        }
    }
    seg <- list(Y = YsegList, Z = ZsegList)
    loc <- list(Y = YlocList, Z = ZlocList)
    
    ## more sanity checks
    stop_if_not(length(YsegList) == dimension)
    nrows <- sapply(seg$Y, nrow); stop_if_not(all(nrows == n))
    nrows <- sapply(seg$Z, nrow); stop_if_not(all(nrows == nbClones))

    if (returnLocus) {
        stop_if_not(length(YlocList) == dimension)
        nrows <- sapply(loc$Y, nrow); stop_if_not(all(nrows == n))
        nrows <- sapply(loc$Z, nrow); stop_if_not(all(nrows == nbClones))
    }
    ## /more sanity checks
    
    if (dimension == 1L) { ## for backward compatibility
        seg$Y <- seg$Y[[1]]
        seg$Z <- seg$Z[[1]]
        if (returnLocus) {
            loc$Y <- loc$Y[[1]]
            loc$Z <- loc$Z[[1]]
        }
    }
        
    res <- list(W = W, segment = seg)
    if (returnLocus) {
        res$locus <- loc
    }
    res
}
