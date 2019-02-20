#' Plot Percentage Variance Explained (PVE) against number of features 
#'
#' @param res Result from [c3co()].
#'
#' @param bestNbLatent Best number of latent profiles.
#'
#' @param ylim A vector that define min and max of y-axis.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @importFrom ggplot2 aes_ ggplot geom_line geom_point geom_vline theme_bw xlab ylim
#' @export
pvePlot <- function(res, bestNbLatent=NULL, ylim=c(0, 1)) {
    pvePlot2(res@config$best, bestNbLatent=bestNbLatent, ylim=ylim)
}


#' Plot Percentage Variance Explained (PVE) against number of features (low-level)
#' 
#' @param dat A data frame, typically the element `best` of the
#'   `config` slot from [c3co()] or [fitC3co()].
#'   
#' @param bestNbLatent Best number of latent profiles.
#'   
#' @param ylim A vector that define min and max of y-axis.
#'   
#' @return A [ggplot2::ggplot] object.
#'   
#' @importFrom ggplot2 aes_ ggplot geom_line geom_point geom_vline
#                      theme_bw xlab ylim
#' @export
pvePlot2 <- function(dat, bestNbLatent=NULL, ylim=c(0, 1)) {
    gg <- ggplot(dat, aes_(x=~nb.feat, y=~PVE))
    gg <- gg + geom_line() + geom_point()
    gg <- gg + ylim(ylim) + xlab("Number of latent profiles")
    gg <- gg + theme_bw()
    if (!is.null(bestNbLatent)) {
        gg <- gg + geom_vline(xintercept=bestNbLatent, lty=2)
    }
    gg
}


#' Plot latent profiles along chromosomes
#' 
#' @param df data.frame Object output from [createZdf()].
#'   
#' @param scalePosToMb A logical, should 'x' positions be scaled to
#'   megabases? Defaults to `FALSE`.
#'   
#' @return A [ggplot2::ggplot] object.
#'   
#' @importFrom ggplot2 ggplot geom_step geom_segment aes facet_grid labeller
#'   label_both theme_bw scale_x_continuous labs scale_y_continuous
#' @export
Zplot <- function(df, scalePosToMb=FALSE) {
    if (scalePosToMb) {
        df$start <- df$start/1e6 
        df$end <- df$end/1e6 
    }
    gg <- ggplot(df, aes_(~start, ~CopyNumber, group=~arch,
                          col=~arch, lty=~arch))
    gg <- gg + geom_step(direction="hv", lwd=1)
    gg <- gg + geom_segment(aes_(xend=~end, yend=~CopyNumber), lwd=1)  ## FIXME: we're actually plotting all the horizontal lines twice
    gg <- gg + facet_grid(stat~chr,
                          labeller=labeller(.cols=label_both))
    gg <- gg + theme_bw()
    gg <- gg + labs(colour = "Subclone", lty="Subclone")
    gg <- gg + scale_x_continuous(name="Genome position")
    gg
}
