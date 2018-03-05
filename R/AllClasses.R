#' Class for the object returned by [positiveFusedLasso()]
#' 
#' @slot Y A list containing the original minor, major and total copy number
#'   signals.
#'   
#' @slot S A list containing the (final) minor, major and total copy number 
#'   estimates for the latent features, as inferred by 
#'   [positiveFusedLasso()].
#'   
#' @slot S0 idem for the initial estimates.
#'   
#' @slot W A matrix containing weights inferred by [positiveFusedLasso()].
#'   
#' @slot mu A numeric containing the intercept terms inferred by
#'   [positiveFusedLasso()].
#'   
#' @slot E A list containing estimates of minor copy number and major copy 
#'   number inferred by [positiveFusedLasso()].
#' 
#' @slot failure A logical value indicating whether the estimation worked or
#'   failed (because of the non-invertibility of \eqn{W^tW}).
#' 
#' @slot params A numeric vector, the input parameters.
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


#' Class for the object create by [fitC3co()]
#'
#' @slot bkp A list of breakpoints for each chromosome.
#'
#' @slot segDat A list that contains segmented data.
#' 
#' @slot config A list containing PVE, BIC and parameters for all tested
#'              models and for the best model.
#'
#' @slot fit A list of [posFused][posFused-class] objects.
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
