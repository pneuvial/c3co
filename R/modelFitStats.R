#' Print statistics of a c3co model
#' 
#' @param Y Original n-by-J matrix of data.
#' 
#' @param Yhat Estimated n-by-J matrix.
#' 
#' @param What Inferred n-by-K weight matrix.
#' 
#' @param Zhatt The transposed version of the inferred n-by-J subclones
#' matrix 'Zhat'.
#' 
#' @return The main statistics of the inferred model.
#'
#' @importFrom matrixStats colDiffs
modelFitStatistics <- function(Y, Yhat, What, Zhatt) {
    stop_if_not(is.matrix(Y), is.matrix(Yhat), is.matrix(What), is.matrix(Zhatt))
    n <- nrow(Y)
    J <- ncol(Y)
    K <- ncol(What)
    stop_if_not(nrow(Yhat) == n, ncol(Yhat) == J)
    stop_if_not(nrow(Zhatt) == J, ncol(Zhatt) == K)
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
    kZ <- sum(colDiffs(Zhatt) != 0)         ## number of breakpoints
    kW <- sum(What != 0)                    ## non null coefs in W
    BIC.Z <- logLik - 1/2*kZ*log(nJ)        ## the one in the manuscript as of July 2017
    BIC.WZ <- logLik - 1/2*(kZ+kW)*log(nJ)  ## the one I (PN) believe we should use
    
    c(BIC=BIC.WZ, PVE=PVE, logLik=logLik, loss=loss)
}

#' @rdname modelFitStats
setMethod("modelFitStats", signature(object = "posFused"), function(object) {
    stats <- modelFitStatistics(Y=object@Y$Y, Yhat=object@E$Y, What=object@W, Zhatt=object@Zt$Z)
    pars <- object@params
    c(pars, stats)
})
