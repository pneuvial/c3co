#'  @title Class for the object return by \code{positiveFusedLasso} function
#'  @slot S A list containing Total copy number, minor copy number and major copy number inferred by \code{positiveFusedLasso}
#'  @slot W A matrix containing weights inferred by \code{positiveFusedLasso}
#'  @slot E A list containing estimates of minor copy number and major copy number inferred by \code{positiveFusedLasso}
#'  @exportClass posFused
#'  
setClass(
  Class = "posFused",
  representation(S = "list", W = "matrix", E = "list")
)
#' @title The method showPosFused
#' @description Print the message that is contained in slots S, W and E.
#' \itemize{
#' \item{showPosFused(this)}{this is an object from class  [\code{\linkS4class{posFused}}].}
#' }
#' @param this an object from the following class: [\code{\linkS4class{posFused}}] 
#' @return nothing
#' @rdname showPosFused
#' @exportMethod showPosFused
setGeneric(
  name = "showPosFused",
  def = function(this) {
    standardGeneric("showPosFused")
  }
)
#' @rdname showPosFused
setMethod(
  f = "showPosFused",
  signature = signature("posFused"),
  def = function(this) {
    cat("Subclones\n")
    cat(utils::str(this@S),"\n")
    cat("Weights\n")
    cat(utils::str(this@W), "\n")
    cat("Estimates\n")
    cat(utils::str(this@E),"\n")
  }
)
#' @title The method Wplot.
#' \itemize{
#' \item{Wplot(this)}{this is an object from class  [\code{\linkS4class{posFused}}].}
#' }
#' @description  Plot the weight matrix
#' @param this an object from the following class: [\code{\linkS4class{posFused}}]
#' @param rownamesW A vector that contains identification of patients
#' @param col A vector that contains colors for the heatmap
#' @param margins A vector margins 
#' @param posLegend position of the legend as for \code{plot}
#' @param listPheno A Matrix that contains details on phenotype for each patient. Could be location or time point of tumors for example
#' @param colsPheno Matrix that containts colors for each type of variable in phenotype
#' @param colLegend colors for clinical data
#' @param labelLegend labels for clinical data
#' @param cexCol size of labels of columns by (default 1.5) 
#' @param ... other paramaters to personalize heatmap (see \code{heatmap.3.R})
#' @return Heatmap of W
#' @rdname Wplot
#' @exportMethod Wplot
setGeneric(
  name = "Wplot",
  def = function(this,rownamesW=NULL, col= NULL, margins=c(5,7), posLegend=NA, listPheno, colsPheno, colLegend, labelLegend, cexCol=1.5,...) {
    standardGeneric("Wplot")
  }
)
#' @rdname Wplot
setMethod(
  f = "Wplot",
  signature = signature("posFused"),
  def = function(this, rownamesW=NULL, col= NULL,margins=c(5,7),posLegend=NA, listPheno, colsPheno, colLegend, labelLegend,cexCol=1.5,...){
    if(is.null(col)){col=grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, 'GnBu'))(100)}
    W <- this@W
    rownames(W) <- rownamesW
    colnames(W) <- sprintf("Subclone %s", letters[1:ncol(W)])
    res.clust = stats::hclust(stats::dist(W),method="ward.D")
    if(!missing(listPheno)){
      if(ncol(listPheno)!=ncol(colsPheno)){
        stop("listPheno and colsPheno must have the same number of columns")
      }
    }
    if(!missing(colsPheno)){
      heatmap.3(W, Rowv=TRUE,dendrogram="row", col=col,scale="none", cexCol=cexCol, cexRow=1.5,margins = margins, key = TRUE,  RowSideColors=t(colsPheno),...)
      
      #heatmap.3(W, Rowv=TRUE, dendrogram="row",  RowSideColors=t(colsPheno), col=col,scale="none", key=TRUE, cexCol=cexCol, cexRow=1.5,margins = c(5,10),...)
      if(!is.na(posLegend)){
        graphics::legend(posLegend,legend=labelLegend, fill=colLegend,border=FALSE, bty="n", y.intersp = 1, cex=1)
      }
    }else{
      heatmap.3(W, Rowv=TRUE,dendrogram="row", col=col,scale="none", cexCol=cexCol, cexRow=1.5,margins = margins, key = TRUE,...)
      #heatmap.3(W, Rowv=TRUE, dendrogram="row", col=col,scale="none", key=TRUE, cexCol=cexCol, cexRow=1.5,margins = margins,...)
    }  
  })



#'  @title Class for the object return by \code{c3co} function
#'  @slot BIC A numeric that is the value of BIC
#'  @slot PVE A numeric that is the Percentage of variation explained
#'  @slot res A \code{posFused} object
#'  @slot param List of parameters (number of subclones, penalty coefficients)
#'  @slot bkp A list of breakpoints for each chromosome
#'  @exportClass c3coClass
#'  
setClass(
  Class = "c3coClass",
  representation(BIC = "numeric", PVE = "numeric", res = "posFused", param = "list", bkp = "list"),
  contains =  "posFused"
)

#' @title The method showC3co
#' @description Print the message that is contained in slots BIC, PVE, res, param and bkp.
#' \itemize{
#' \item{showC3coClass(this)}{this is an object from class  [\code{\linkS4class{c3coClass}}].}
#' }
#' @param this an object from the following class: [\code{\linkS4class{c3coClass}}] 
#' @return nothing
#' @rdname showC3coClass
#' @exportMethod showC3coClass
setGeneric(
  name = "showC3coClass",
  def = function(this) {
    standardGeneric("showC3coClass")
  }
)
#' @rdname showC3coClass
setMethod(
  f = "showC3coClass",
  signature = signature("c3coClass"),
  def = function(this) {
    cat("BIC\n")
    cat(this@BIC,"\n")
    cat("PVE\n")
    cat(this@PVE, "\n")
    cat("posFused\n")
    cat(utils::str(this@res),"\n")
    cat("param\n")
    cat(utils::str(this@param),"\n")
    cat("breakpoints\n")
    cat(utils::str(this@bkp),"\n")
  }
)


#' @title The method createZdf
#' @description Create a data frame to plot Subclones
#' \itemize{
#' \item{createZdf(this)}{this is an object from class  [\code{\linkS4class{c3coClass}}].}
#' }
#' @title The method createZdf
#' @param this an object from the following class: [\code{\linkS4class{c3coClass}}] 
#' @param minMaxPos Matrix that contains min and max position for each chromosome
#' @param chromosomes A vector that contains the focused chromosomes
#' @param var TCN, Minor or Major 
#' @return A data frame to plot Latent profiles with ggplot
#' @rdname createZdf
#' @exportMethod createZdf
setGeneric(
  name = "createZdf",
  def = function(this, minMaxPos, chromosomes, var="TCN") {
    standardGeneric("createZdf")
  }
)
#' @rdname createZdf
setMethod(
  f = "createZdf",
  signature = signature("c3coClass"),
  def = function(this, minMaxPos, chromosomes, var="TCN"){
    lengthCHR <- sapply(this@bkp, length)
    start <- c(1,cumsum(lengthCHR)+1)
    
    df.CHR <- do.call(rbind, lapply(chromosomes, function(cc){
      bb <- c(minMaxPos[cc,"minPos"], this@bkp[[cc]],minMaxPos[cc,"maxPos"])
      bb <- as.numeric(bb)
      if(var=="TCN"){
        zz <- rbind(this@res@S$Z[start[cc],],this@res@S$Z[start[cc]:(start[cc+1]-1),],this@res@S$Z[start[cc+1]-1,])
        
      }else if(var=="Minor"){
        zz <- rbind(this@res@S$Z1[start[cc],],this@res@S$Z1[start[cc]:(start[cc+1]-1),],this@res@S$Z1[start[cc+1]-1,])
        
      }else if(var=="Major"){
        zz <- rbind(this@res@S$Z2[start[cc],],this@res@S$Z2[start[cc]:(start[cc+1]-1),],this@res@S$Z2[start[cc+1]-1,])
      }else{
        stop("var must be TCN, Minor or Major")
      }
      
      arch <- factor(rep(letters[1:ncol(zz)], each=length(bb))) 
      zz <- c(zz)
      data.frame(position=bb,CN=zz, arch=arch, chr=factor(sprintf("chromosome %s",cc))) 
    }))
    return(df.CHR)
  })



