#' Print statistics of a c3co model
#' 
#' @param Y Original n-by-J matrix of data.
#' 
#' @param Yhat Estimated n-by-J matrix.
#' 
#' @param What Inferred n-by-K weight matrix.
#' 
#' @param Zhat Inferred n-by-J subclones matrix.
#' 
#' @return The main statistics of the inferred model.
#'
#' @importFrom matrixStats colDiffs
modelFitStatistics <- function(Y, Yhat, What, Zhat) {
    stop_if_not(is.matrix(Y), is.matrix(Yhat), is.matrix(What), is.matrix(Zhat))
    n <- nrow(Y)
    J <- ncol(Y)
    K <- ncol(What)
    stop_if_not(nrow(Yhat) == n, ncol(Yhat) == J)
    stop_if_not(nrow(Zhat) == J, ncol(Zhat) == K)
    stop_if_not(nrow(What) == n, ncol(What) == K)

    nJ <- n*J
    
    ## PVE
    totSS <- sum(sweep(Y, MARGIN = 1L, STATS = rowMeans(Y), FUN = `-`)^2)
    resSS <- sum((Y - Yhat)^2)
    PVE <- 1 - resSS / totSS

    loss <- resSS/nJ
    
    ## log-likelihood
    logLik <- -(nJ*log(loss) + nJ*(1+log(2*pi)))/2  ## the last bit is indep of model complexity
    ## why *1+*?
    
    ## BIC penalty    
    kZ <- sum(colDiffs(Zhat) != 0)          ## number of breakpoints
    kW <- sum(What != 0)                    ## non null coefs in W
    BIC.Z <- logLik - 1/2*kZ*log(nJ)        ## the one in the manuscript as of July 2017
    BIC.WZ <- logLik - 1/2*(kZ+kW)*log(nJ)  ## the one I (PN) believe we should use
    
    c(BIC=BIC.WZ, PVE=PVE, logLik=logLik, loss=loss)
}

#' @rdname modelFitStats
setMethod("modelFitStats", signature(object = "posFused"), function(object) {
    stats <- modelFitStatistics(object@Y$Y, object@E$Y, object@W, object@S$Z)
    pars <- object@params
    c(pars, stats)
})
