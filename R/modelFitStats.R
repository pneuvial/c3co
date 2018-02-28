#' Print statistics of a c3co model
#' 
#' @param Y Original matrix of data
#' @param Yhat Estimates 
#' @param What Inferred weight matrix
#' @param Zhat Inferred subclones matrix
#' @return The main statistics of the inferred model.
modelFitStatistics <- function(Y, Yhat, What, Zhat) {
    n <- nrow(Y)
    nseg <- ncol(Y)
    nnS <- n*nseg
    
    ## PVE
    totSS <- sum(sweep(Y, MARGIN = 1L, STATS = rowMeans(Y), FUN = `-`)^2)
    resSS <- sum((Y - Yhat)^2)
    PVE <- 1 - resSS / totSS

    loss <- resSS/nnS
    
    ## log-likelihood
    logLik <- -(nnS*log(loss) + nnS*(1+log(2*pi)))/2  ## the last bit is indep of model complexity
    ## why *1+*?
    
    ## BIC penalty    
    kZ <- sum(apply(Zhat, MARGIN=2L, FUN=diff) != 0)  ## number of breakpoints
    kW <- sum(What != 0)                              ## non null coefs in W
    BIC.Z <-  logLik - 1/2*kZ*log(nnS)                ## the one in the manuscript as of July 2017
    BIC.WZ <-  logLik - 1/2*(kZ+kW)*log(nnS)          ## the one I (PN) believe we should use
    
    c(BIC=BIC.WZ, PVE=PVE, logLik=logLik, loss=loss)
}

#' @rdname modelFitStats
setMethod("modelFitStats", signature(object = "posFused"), function(object) {
    stats <- modelFitStatistics(object@Y$Y, object@E$Y, object@W, object@S$Z)
    pars <- object@params
    c(pars, stats)
})
