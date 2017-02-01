#'  @title Class for the object return by \code{positiveFusedLasso} function
#'  @slot S A list containing Total copy number, minor copy number and major copy number inferred by \code{positiveFusedLasso}
#'  @slot W A matrix containing weights inferred by \code{positiveFusedLasso}
#'  @slot E A list containing estimates of minor copy number and major copy number inferred by \code{positiveFusedLasso}
#'  @slot BIC A numeric containing BIC value
#'  @slot PVE A numeric containing PVE value (Percentage of variation explained)
#'  @slot param A list of parameters (number of subclones, penalty coefficients)
#'  @exportClass posFused
#'  
setClass(
  Class = "posFused",
  representation(S = "list", W = "matrix", E = "list", BIC = "numeric", PVE = "numeric", param = "list")
)

#'  @title Class for the object create by \code{c3coFit} function
#'  @slot bkp A list of breakpoints for each chromosome
#'  @slot segDat A list that contains segmented data
#'  @slot fit A List of [\code{\linkS4class{posFused}}] objects
#'  @exportClass c3coClass
#'
setClass("c3coFit", representation(bkp = "list", segDat = "list", fit = "list"))
