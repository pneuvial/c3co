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
        cat("Intercept vector:\n")
        str(object@mu)
        cat("Original signal:")
        str(object@Y)
        cat("Copy-number estimates:\n")
        str(object@E)
    }
)

#' Compute statistics of a c3co model
#'
#' @param object An object of class [posFused][posFused-class].
#' @rdname modelFitStats
#' @exportMethod modelFitStats
setGeneric("modelFitStats", function(object) {
    standardGeneric("modelFitStats")
})

#' Plot the weight matrix
#'
#' @param this An object of class [c3coFit][c3coFit-class].
#'
#' @param idxBest A integer, the number of latent features.
#'
#' @param rownamesW A vector that contains identification of patients.
#'
#' @param col A vector that contains colors for the heatmap.
#'
#' @param margins A vector margins.
#'
#' @param posLegend Position of the legend to be passed to [graphics::plot()].
#'
#' @param listPheno A matrix that contains details on phenotype for each
#'        patient. Could be location or time point of tumors for example.
#'
#' @param colsPheno Matrix that contains colors for each type of variable
#'        in phenotype.
#'
#' @param colLegend Colors for clinical data.
#'
#' @param labelLegend Labels for clinical data.
#'
#' @param cexCol Size of labels of columns by (default 1.5).
#'
#' @param ... Other parameters to personalize heatmap (see [heatmap.3()]).
#'
#' @return A heatmap of `W`.
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
        if (is.null(col)) {
           col <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
        }
        W <- this@fit[[idxBest]]@W
        rownames(W) <- rownamesW
        colnames(W) <- sprintf("Subclone %s", letters[1:ncol(W)])
        if (!missing(listPheno)) {
            if (ncol(listPheno) != ncol(colsPheno)) {
                stop("'listPheno' and 'colsPheno' must have the same number of columns")
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
        cat("Configuration:\n")
        str(object@config)
        cat("Results of positive fused lasso:\n")
        cat("List of", length(object@fit), "objects of class 'posFused'")
    }
)

#' Create a data frame to plot subclones
#'
#' @param this An object of class [c3coFit][c3coFit-class].

#' @param chromosomes A vector that contains the focused chromosomes.

#' @param var `"TCN"`, `"Minor"` or `"Major"`.

#' @param idxBest A integer that is the best fitting of the data.

#' @param verbose A logical value indicating whether to print extra information.
#'   Defaults to `FALSE`.

#' @return A data frame to plot Latent profiles with \pkg{ggplot2}.

#' @rdname createZdf
#' @exportMethod createZdf
setGeneric(
    name = "createZdf",
    def = function(this, chromosomes, idxBest, var="TCN", verbose=FALSE) {
        standardGeneric("createZdf")
    }
)
#' @rdname createZdf
setMethod(
    f = "createZdf",
    signature = signature("c3coFit"),
    definition = function(this, chromosomes, idxBest, var=c("TCN", "Minor", "Major"), verbose=FALSE) {
        var <- match.arg(var, several.ok=TRUE)
        labs <- list(TCN="Z", Minor="Z1", Major="Z2")
        bkp <- this@bkp
        lengthCHR <- sapply(bkp, FUN=function(x) length(x)-1) ## '-1' because 'bkp' includes first and last position on chr
        idx <- c(1, cumsum(lengthCHR) + 1)
        fitZ <- this@fit[[idxBest]]@S
        nbarch <- ncol(fitZ$Z)

        dfList <- list()
        configs <- expand.grid(var=var, chr=chromosomes)
        if (verbose) mprint(configs)
        for (kk in 1:nrow(configs)) {
            vv <- as.character(configs[kk, "var"])
            cc <- configs[kk, "chr"]
            lab <- labs[[vv]]
            Z <- fitZ[[lab]]
            if (is.null(Z)) next
            idxCC <- idx[cc]:(idx[cc+1]-1)  ## indices of observations in current chr
            Z <- Z[idxCC,, drop=FALSE]
##            Z <- rbind(Z, Z[nrow(Z), ])  ## ad hoc: repeat last segment so that geom_step plots it?
            dim(Z) <- NULL ## convert matrix to vector
            bb <- bkp[[cc]]
            bb <- as.numeric(bb)
            nbseg <- length(bb) - 1 ## '-1' because 'bb' includes first and last position on chr 
            start <- bb[1:nbseg]
            end <- bb[1:nbseg+1]
            arch <- factor(rep(letters[1:nbarch], each=nbseg))
            datCC <- data.frame(chr=cc, start=start, arch=arch, end=end, CopyNumber=Z, stat=vv)
            dfList[[kk]] <- datCC
        }
        df.CHR <- do.call(rbind, args=dfList)
        df.CHR
    }
)
