#' Class for the object returned by \code{positiveFusedLasso} function
#'
#' @slot S A list containing the (final) minor, major and total copy number
#'       estimates for the latent features, as inferred by
#'       \code{\link{positiveFusedLasso}}.
#'
#' @slot S0 idem for the initial estimates.
#'
#' @slot W A matrix containing weights inferred by \code{positiveFusedLasso}.
#'
#' @slot mu A numeric containing the intercept terms inferred by \code{positiveFusedLasso}.
#'
#' @slot E A list containing estimates of minor copy number and major copy
#'       number inferred by \code{positiveFusedLasso}.
#'       
#' @exportClass posFused
setClass(
    Class = "posFused",
    representation = representation(
        Y  = "list",        ## original signal
        S  = "list",        ## subclones (aka Z)
        S0 = "list",        ## initial estimates of subclones (aka Z0)
        W  = "matrix",      ## weights
        mu = "numeric",     ## intercept
        E  = "list",        ## signal reconstruction (aka Yhat)
        params  = "numeric", 
        failure = "logical"
    )
)


#' Class for the object create by \code{c3coFit} function
#'
#' @slot bkp A list of breakpoints for each chromosome.
#'
#' @slot segDat A list that contains segmented data.
#' 
#' @slot config A list containing PVE, BIC and parameters for all tested models 
#'                and for the best model
#'
#' @slot fit A List of [\code{\linkS4class{posFused}}] objects.
#'
#' @exportClass c3coFit
setClass(
    Class = "c3coFit",
    representation = representation(
        bkp = "list",
        segDat = "list",
        config = "list",
        fit = "list"
    )
)
