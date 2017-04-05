#' Function to plot PVE
#'
#' @param res result from [\code{posFused}]   
#' @param bestNbLatent best number of latent profiles.
#' @param ylim a vector that define min and max of y-axis
#' @return PVE curve
#'
#' @importFrom ggplot2 aes_ ggplot geom_line geom_point geom_vline theme_bw xlab ylim
#' @export
pvePlot <- function(res, bestNbLatent=NULL, ylim=c(0, 1)) {
    PVEs <- sapply(res, FUN = function(rr) rr@PVE)
    nb.arch <- sapply(res, FUN = function(rr) rr@param$nb.arch)
    df.PVE <- data.frame(PVE=PVEs, nb.arch=nb.arch)
    gg <- ggplot(df.PVE, aes_(x=~nb.arch, y=~PVE)) + geom_line() + geom_point() + theme_bw() + ylim(ylim) + xlab("Number of latent profiles")
    if (!is.null(bestNbLatent)) {
        gg <- gg + geom_vline(xintercept=bestNbLatent, lty=2)
    }
    gg  
}



#' Plot latent profiles along chromosomes
#'
#' @export
#' @param df data.frame object output from \code{createZdf}
#' @return plot latent profiles along chromosomes
#' @importFrom ggplot2 ggplot geom_step aes_ facet_grid labeller label_both theme_bw scale_x_continuous labs
Zplot <- function(df) {
    gg <- ggplot(df)
    gg <- gg + geom_step(aes_(~position, ~CopyNumber, group=~arch, col=~arch, lty=~arch), direction="hv", lwd=1)
    gg <- gg + facet_grid(stat~chr, scales="free", labeller=labeller(.cols=label_both))
    gg <- gg + theme_bw()
    gg <- gg + labs(colour = "Subclone", lty="Subclone")
    gg <- gg + scale_x_continuous(name="Genome position")
    gg
}
