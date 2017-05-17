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
#' @slot E A list containing estimates of minor copy number and major copy
#'       number inferred by \code{positiveFusedLasso}.
#'       
#' @exportClass posFused
setClass(
  Class = "posFused",
  representation = representation(
    S = "list",
    S0 = "list",
    W = "matrix",
    E = "list"
  )
)


#' Class for the object create by \code{c3coFit} function
#'
#' @slot bkp A list of breakpoints for each chromosome.
#'
#' @slot segDat A list that contains segmented data.
#' 
#' @slot config A data.frame containing PVE, BIC and parameters of models
#'
#' @slot fit A List of [\code{\linkS4class{posFused}}] objects.
#'
#' @exportClass c3coFit
setClass(
  Class = "c3coFit",
  representation = representation(
    bkp = "list",
    segDat = "list",
    config = "data.frame",
    fit = "list"
  )
)
