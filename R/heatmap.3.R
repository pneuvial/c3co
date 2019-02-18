## Source: Benjamin Sadacca, Institut Curie
## ## trim.heatmap
#' @importFrom stats quantile
trim.heatmap <- function(data, trim) {
    ## data <- data - mean(data, na.rm = TRUE)
    data <- t(scale(t(data)))
    q <- quantile(data, probs = c((1 - trim), trim), na.rm = TRUE)
    data[data < q[1]] <- q[1]
    data[data > q[2]] <- q[2]
    maxi <- max(data, na.rm = TRUE)
    mini <- min(data, na.rm = TRUE)
    data[!is.na(data) & data > 0] <-  data[!is.na(data) &  data > 0]/maxi
    data[!is.na(data) & data < 0] <- -data[!is.na(data) &  data < 0]/mini
    data
}

#' Heatmap function
#'
#' @param x, matrix
#' @param Rowv = TRUE,
#' @param Colv = if (symm) "Rowv" else TRUE,
#' @param distfun = dist,
#' @param hclustfun = hclust,
#' @param dendrogram = c("both","row", "column", "none"),
#' @param symm = FALSE,
#' @param scale = c("none","row", "column"),
#' @param na.rm = TRUE,
#' @param revC = identical(Colv,"Rowv"),
#' @param add.expr, expression that will be evaluated after the call to image. Can be used to add components to the plot.
#' @param breaks, numeric, either a numeric vector indicating the splitting points for binning x into colors, or a integer number of break points to be used, in which case the break points will be spaced equally between range(x). DEFAULT: 16 when not specified.
#' @param symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
#' @param col = "heat.colors",
#' @param colsep, color for columns of separation.
#' @param rowsep, color for lines of separation.
#' @param sepcolor = "white",
#' @param sepwidth = c(0.05, 0.05),
#' @param cellnote, optional) matrix of character strings which will be placed within each color cell, e.g. cell labels or p-value symbols.
#' @param notecex = 1,
#' @param notecol = "cyan",
#' @param na.color = par("bg"),
#' @param trace = c("none", "column","row", "both"),
#' @param tracecol = "cyan",
#' @param hline = median(breaks),
#' @param vline = median(breaks),
#' @param linecol = tracecol,
#' @param margins = c(5,5),
#' @param ColSideColors, (optional) character vector of length ncol(x) containing the color names for a horizontal side bar that may be used to annotate the columns of x.
#' @param RowSideColors, (optional) character vector of length nrow(x) containing the color names for a horizontal side bar that may be used to annotate the rows of x.
#' @param side.height.fraction = 0.3, the width-to-height ratio of the heatmap.
#' @param cexRow = 0.2 + 1/log10(nr),
#' @param cexCol = 0.2 + 1/log10(nc),
#' @param cexMain = 1.5, size of main
#' @param cexKey = 1, size of name of key
#' @param cexColorKey = 1, ## size of key
#' @param labRow = NULL,
#' @param labCol = NULL,
#' @param key = TRUE,
#' @param keysize = 1.5,
#' @param density.info = c("none", "histogram", "density"),
#' @param denscol = tracecol,
#' @param symkey = max(x < 0, na.rm = TRUE) || symbreaks,
#' @param densadj = 0.25,
#' @param main = NULL,
#' @param xlab = NULL,
#' @param ylab = NULL,
#' @param lmat = NULL,
#' @param lhei = NULL,
#' @param lwid = NULL,
#' @param ColSideColorsSize = 1,
#' @param RowSideColorsSize = 1,
#' @param KeyValueName = "", name of key
#' @param ... additional parameters of [graphics::image()].
#' @return A Heatmap
#'
#' @importFrom graphics abline axis hist image layout lines mtext par plot plot.new rect text title
#' @importFrom stats as.dendrogram density dist hclust median order.dendrogram reorder
#' @importFrom matrixStats colSds rowSds
#' @export
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both", "row", "column", "none"),
                      symm = FALSE,
                      scale = c("none", "row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv, "Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column", "row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5, 5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      cexMain = 1.5, ## add me
                      cexKey=1, ## add me
                      cexColorKey=1, ## add me
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="", ...) {

    invalid <- function(x) {
        if (missing(x) || is.null(x) || length(x) == 0)
            return(TRUE)
        if (is.list(x))
            return(all(sapply(x, FUN=invalid)))
        else if (is.vector(x))
            return(all(is.na(x)))
        else return(FALSE)
    }

    x <- as.matrix(x)
    scale01 <- function(x, low = 0, high = 1) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
                "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || (!inherits(Rowv, "dendrogram") && is.na(Rowv)))
        Rowv <- FALSE
    if (is.null(Colv) || (!inherits(Colv, "dendrogram") && is.na(Colv)))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("'x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("'margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                     c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                    dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                     c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: 'Colv' is FALSE, while dendrogram is `",
                    dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }


    ddc <- Colv
    colInd <- 1:nc
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, ]  ##FIXME: drop=FALSE or drop=TRUE?
    x.unscaled <- x
    cellnote <- cellnote[rowInd, ]  ##FIXME: drop=FALSE or drop=TRUE?
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, MARGIN=1L, STATS=rm)
        retval$rowSDs <- sx <- rowSds(x, na.rm = na.rm)
        x <- sweep(x, MARGIN=1L, STATS=sx, FUN=`/`)
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, MARGIN=2L, STATS=rm)
        retval$colSDs <- sx <- colSds(x, na.rm = na.rm)
        x <- sweep(x, MARGIN=2L, STATS=sx, FUN=`/`)
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(from=0, to=1,
                          length.out = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(from=-extreme, to=extreme, length.out = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (is.function(col))
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)

        if (!missing(ColSideColors)) {
            #if (!is.matrix(ColSideColors))
            #stop("'ColSideColors' must be a matrix")
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA_real_, 1), lmat[2, ] + 1)
            #lhei <- c(lhei[1], 0.2, lhei[2])
            lhei <- c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
        }

        if (!missing(RowSideColors)) {
            #if (!is.matrix(RowSideColors))
            #stop("'RowSideColors' must be a matrix")
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA_real_, times=nrow(lmat) - 1), 1), lmat[, 2] + 1)
            #lwid <- c(lwid[1], 0.2, lwid[2])
            lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }

    if (length(lhei) != nrow(lmat))
        stop("'lhei' must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("'lwid' must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)) {
            par(mar = c(margins[1], 0, 0, 0.5))
            image(x=rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc <- t(RowSideColors[, rowInd, drop=FALSE])
            rsc.colors <- matrix()  ## FIXME: == matrix(NA) - inefficient data type
            rsc.names <- names(table(rsc))
            rsc.i <- 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] <- rsc.name
                rsc[rsc == rsc.name] <- rsc.i
                rsc.i <- rsc.i + 1
            }
            rsc <- matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(x=t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(rownames(RowSideColors)) > 0) {
                axis(side=1L, at=0:(dim(rsc)[2] - 1)/max(1, (dim(rsc)[2] - 1)), labels=rownames(RowSideColors), las = 2, tick = FALSE)
            }
        }
    }

    if (!missing(ColSideColors)) {

        if (!is.matrix(ColSideColors)) {
            par(mar = c(0.5, 0, 0, margins[2]))
            image(x=cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc <- ColSideColors[colInd, , drop=FALSE]
            csc.colors <- matrix()  ## FIXME: == matrix(NA) - inefficient data type
            csc.names <- names(table(csc))
            csc.i <- 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] <- csc.name
                csc[csc == csc.name] <- csc.i
                csc.i <- csc.i + 1
            }
            csc <- matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(side=2L, at=0:(dim(csc)[2] - 1)/max(1, (dim(csc)[2] - 1)), labels=colnames(ColSideColors), las = 2, tick = FALSE)
            }
        }
    }

    par(mar = c(margins[1], 0, 0.5, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]  ## FIXME: drop=FALSE or drop=TRUE?
    }
    else iy <- 1:nr
    image(x=1:nc, y=1:nr, z=x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) && anyNA(x)) { # load library(gplots)
        mmat <- ifelse(is.na(x), 1, NA_real_)
        image(x=1:nc, y=1:nr, z=mmat, axes = FALSE, xlab = "", ylab = "",
                        col = na.color, add = TRUE)
    }
    axis(side=1L, at=1:nc, labels = labCol, las = 1, line = -0.5, tick = 0,
                   cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(side=4L, at=iy, labels = labRow, las = 2, line = -0.5, tick = 0,
                   cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, times=length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, times=csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                                 lty = 2)
            }
            xv <- rep(i, times=nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, times=ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
                       col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i")
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = cexMain * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- 0
            max.raw <- 1
        }

        z <- seq(from=min.raw, to=max.raw, length.out = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
                        xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(side=1L, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2, cex=cexKey)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                            lwd = 1)
            axis(side=2L, at = pretty(dens$y)/max(dens$y) * 0.95, labels=pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                            col = denscol)
            axis(side=2L, at = pretty(hy)/max(hy) * 0.95, labels=pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("percentage of
clone in sample", cex.main=cexColorKey)
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                    high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}
