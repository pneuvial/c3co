#' Print statistics of a c3co model
#' 
#' @param Y Original n-by-J matrix of data for n samples and J segments.
#' 
#' @param Yhat Estimated n-by-J matrix for n samples and J segments.
#' 
#' @param What Inferred n-by-K weight matrix for n samples and K subclones.
#' 
#' @param Zhatt The transposed version of the inferred n-by-J subclones
#' matrix 'Zhat' for n samples and J segments.
#' 
#' @return A named numeric vector with the main statistics of the inferred model:
#'  * `BIC`   : Bayesian Information Criterion
#'  * `PVE`   : Percentage of Variation Explained (Nowak et al., 2011)
#'  * `logLik`: Log Likelihood
#'  * `loss`  : Loss, which equals \eqn{RSS / (n * J)},
#'              where RSS is the residual sum of squares
#'
#' @references Nowak, G., Hastie, T., Pollack, J. R., & Tibshirani, R. (2011).
#'             A fused lasso latent feature model for analyzing multi-sample
#'             aCGH data. Biostatistics, 12(4), 776-791
#'
#' @importFrom matrixStats colDiffs
modelFitStatistics <- function(Y, Yhat, What, Zhatt) {
    stop_if_not(is.matrix(Y), is.matrix(Yhat), is.matrix(What), is.matrix(Zhatt))
    n <- nrow(Y)     ## Number of samples
    J <- ncol(Y)     ## Number of segments
    K <- ncol(What)  ## Number of subclones
    stop_if_not(nrow(Yhat) == n, ncol(Yhat) == J)
    stop_if_not(nrow(Zhatt) == J, ncol(Zhatt) == K)
    stop_if_not(nrow(What) == n, ncol(What) == K)

    nJ <- n*J
    
    ## Percentage of Variation Explained (PVE)
    ## Comment: We have observed that PVE is sometimes negative [#19].
    ##          This can only happen if below resSS > totSS.  So, when
    ##          does that occur? /HB 2019-02-21
    totSS <- sum((Y - rowMeans(Y))^2)
    resSS <- sum((Y - Yhat)^2)
    PVE <- 1 - resSS / totSS

    ## FIXME: Clarified the calculation of the total sum of squares to
    ##        rule out a thinko in the usage of sweep().  Keeping both
    ##        for now and using run-time checks assert they are the same.
    ##        /HB 2019-02-21
    totSS_old <- sum(sweep(Y, MARGIN = 1L, STATS = rowMeans(Y), FUN = `-`)^2)
    stop_if_not(totSS == totSS_old)

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
