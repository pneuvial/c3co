#' Class for the object returned by [positiveFusedLasso()]
#' 
#' @slot Y A list containing the original minor, major and total copy number
#'   signals.
#'   
#' @slot Zt A list containing the (final) minor, major and total copy number 
#'   estimates for the latent features, as inferred by 
#'   [positiveFusedLasso()].
#'   
#' @slot Z0t idem for the initial estimates.
#'   
#' @slot W A matrix containing weights inferred by [positiveFusedLasso()].
#'   
#' @slot mu A numeric containing the intercept terms inferred by
#'   [positiveFusedLasso()].
#'   
#' @slot E A list containing estimates of minor copy number and major copy 
#'   number inferred by [positiveFusedLasso()].
#' 
#' @slot params A numeric vector, the input parameters.
#' 
#' @slot failure A logical value indicating whether the estimation worked or
#'   failed (because of the non-invertibility of \eqn{W^tW}).
#' 
#' @slot converged,iterations Did the model fit converge and how many
#'   iterations did it run?
#'   
#' @exportClass posFused
setClass(
    Class = "posFused",
    representation = representation(
        Y   = "list",           ## original signal
        Zt  = "list",           ## subclones (was S)
        Z0t = "list",           ## initial estimates of subclones (was S0)
        W   = "matrix",         ## weights
        mu  = "numeric",        ## intercept
        E   = "list",           ## signal reconstruction (aka Yhat)
        params  = "numeric",
        failure = "logical",
	converged = "logical",  ## Did the model fit converge?
	iterations = "integer"  ## How many iterations did it run?
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
