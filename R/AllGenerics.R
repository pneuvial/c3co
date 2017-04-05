#' @importFrom utils str
setMethod(
    f = "show",
    signature = signature("posFused"),
    definition = function(object) {
        cat("Object of class 'posFused':\n")
        cat("Subclones:")
        str(object@S)
        cat("Weight matrix:\n")
        str(object@W)
        cat("Copy-number estimates:\n")
        str(object@E)
        cat("Bayesian Information Criterion (BIC):", object@BIC, "\n")
        cat("Percentage of Variance Explained (PVE):", object@PVE, "\n")
        cat("Parameters:\n")
        str(object@param)
    }
)

#' Plot the weight matrix
#'
#' @param this An object of class \code{\linkS4class{c3coFit}}
#'
#' @param idxBest A integer, the number of latent features
#'
#' @param rownamesW A vector that contains identification of patients
#'
#' @param col A vector that contains colors for the heatmap
#'
#' @param margins A vector margins
#'
#' @param posLegend position of the legend to be passed to \code{plot}
#'
#' @param listPheno A matrix that contains details on phenotype for each
#'        patient. Could be location or time point of tumors for example
#'
#' @param colsPheno Matrix that containts colors for each type of variable
#'        in phenotype
#'
#' @param colLegend colors for clinical data
#'
#' @param labelLegend labels for clinical data
#'
#' @param cexCol size of labels of columns by (default 1.5)
#'
#' @param ... other paramaters to personalize heatmap (see \code{heatmap.3})
#'
#' @return a heatmap of W
#'
#' @rdname Wplot
#' @exportMethod Wplot
setGeneric(
    name = "Wplot",
    def = function(this, idxBest, rownamesW=NULL, col=NULL, margins=c(5, 7), posLegend=NA, listPheno, colsPheno, colLegend, labelLegend, cexCol=1.5, ...) {
        standardGeneric("Wplot")
    }
)
#' @importFrom graphics legend
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @rdname Wplot
setMethod(
    f = "Wplot",
    signature = signature("c3coFit"),
    definition = function(this, idxBest, rownamesW=NULL, col=NULL, margins=c(5, 7), posLegend=NA, listPheno, colsPheno, colLegend, labelLegend, cexCol=1.5, ...) {
        if (is.null(col)) { col <- colorRampPalette(brewer.pal(9, 'GnBu'))(100) }
        W <- this@fit[[idxBest]]@W
        rownames(W) <- rownamesW
        colnames(W) <- sprintf("Subclone %s", letters[1:ncol(W)])
        if (!missing(listPheno)) {
            if (ncol(listPheno) != ncol(colsPheno)) {
                stop("listPheno and colsPheno must have the same number of columns")
            }
        }
        if (!missing(colsPheno)) {
            heatmap.3(W, Rowv=TRUE, dendrogram="row", col=col, scale="none", cexCol=cexCol, cexRow=1.5, margins = margins, key = TRUE, RowSideColors=t(colsPheno), ...)
            
            #heatmap.3(W, Rowv=TRUE, dendrogram="row", RowSideColors=t(colsPheno), col=col, scale="none", key=TRUE, cexCol=cexCol, cexRow=1.5, margins = c(5, 10), ...)
            if (!is.na(posLegend)) {
                legend(posLegend, legend=labelLegend, fill=colLegend, border=FALSE, bty="n", y.intersp = 1, cex=1)
            }
        } else {
            heatmap.3(W, Rowv=TRUE, dendrogram="row", col=col, scale="none", cexCol=cexCol, cexRow=1.5, margins = margins, key = TRUE, ...)
            #heatmap.3(W, Rowv=TRUE, dendrogram="row", col=col, scale="none", key=TRUE, cexCol=cexCol, cexRow=1.5, margins = margins, ...)
        }
    })


#' @importFrom utils str
setMethod(
    f = "show",
    signature = signature("c3coFit"),
    definition = function(object) {
        cat("Object of class 'c3coFit':\n")
        cat("Breakpoints:\n")
        str(object@bkp)
        cat("Segmented data:\n")
        str(object@segDat)
        cat("Results of positive fused lasso:\n")
        cat("List of", length(object@fit), "objects of class 'posFused'")
    }
)

#' Create a data frame to plot Subclones
#'
#' @param this An object of class \code{\linkS4class{c3coFit}}
#' @param minMaxPos Matrix containing min and max position for each chromosome
#' @param chromosomes A vector that contains the focused chromosomes
#' @param var TCN, Minor or Major
#' @param idxBest a integer that is the best fitting of the data
#' @return A data frame to plot Latent profiles with ggplot
#' @rdname createZdf
#' @exportMethod createZdf
setGeneric(
    name = "createZdf",
    def = function(this, minMaxPos, chromosomes, var="TCN", idxBest) {
        standardGeneric("createZdf")
    }
)
#' @rdname createZdf
setMethod(
    f = "createZdf",
    signature = signature("c3coFit"),
    definition = function(this, minMaxPos, chromosomes, var=c("TCN", "Minor", "Major"), idxBest) {
        var <- match.arg(var, several.ok=TRUE)
        labs <- list("TCN"="Z", "Minor"="Z1", "Major"="Z2")
        bkp <- this@bkp
        lengthCHR <- sapply(bkp, FUN=length)
        start <- c(1, cumsum(lengthCHR)+1)
        fitZ <- this@fit[[idxBest]]@S
        
        dfList <- list()
        configs <- expand.grid(var=var, chr=chromosomes)
        print(configs)
        for (kk in 1:nrow(configs)) {
            vv <- configs[kk, "var"]
            cc <- configs[kk, "chr"]
            lab <- labs[[vv]]
            Z <- fitZ[[lab]]
            
            bb <- c(minMaxPos[cc, "minPos"], bkp[[cc]], minMaxPos[cc, "maxPos"])
            bb <- as.numeric(bb)
            zz <- rbind(Z[start[cc], ],
                        Z[start[cc]:(start[cc+1]-1), ],
                        Z[start[cc+1]-1, ])
            arch <- factor(rep(letters[1:ncol(zz)], each=length(bb)))
            zz <- c(zz)
            datCC <- data.frame(position=bb, CopyNumber=zz, arch=arch,
                                chr=cc, stat=vv)
            dfList[[kk]] <- datCC
        }
        df.CHR <- do.call(rbind, args=dfList)
        return(df.CHR)
    }
)
