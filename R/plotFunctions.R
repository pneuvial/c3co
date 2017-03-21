#' Function to plot PVE
#'
#' @export
#' @param res result from [\code{posFused}]   
#' @param bestNbLatent best number of latent profiles.
#' @param ylim a vector that define min and max of y-axis
#' @return PVE curve
pvePlot <- function(res, bestNbLatent=NULL, ylim=c(0,1)){
  PVEs <- sapply(res, function (rr) rr@PVE)
  nb.arch <- sapply(res, function (rr) rr@param$nb.arch)
  df.PVE <- data.frame(PVE=PVEs, nb.arch=nb.arch)
  gg <- ggplot2::ggplot(df.PVE,ggplot2::aes_(x=~nb.arch, y=~PVE))+ ggplot2::geom_line()+ggplot2::geom_point()+ggplot2::theme_bw()+ggplot2::ylim(ylim)+ggplot2::xlab("Number of latent profiles")
  if(!is.null(bestNbLatent)){
    gg <- gg+ggplot2::geom_vline(xintercept=bestNbLatent, lty=2)
  }
  gg  
}



#' Function to plot Latent profiles
#'
#' @export
#' @param df data.frame object output from \code{createZdf}
#' @return plot of Latent profiles
Zplot <- function(df, ylab) {
  gg <- ggplot2::ggplot(df)
  gg <- gg+ggplot2::geom_step(ggplot2::aes_(~position, ~CopyNumber, group=~arch, col=~arch,lty=~arch), direction="hv", lwd=1)
  gg <- gg+ggplot2::facet_grid(stat~chr, scales="free", labeller=labeller(.cols=label_both))
  gg <- gg+ggplot2::theme_bw()
  gg <- gg+ggplot2::labs(colour = "Subclone", lty="Subclone")
  gg <- gg+ggplot2::scale_x_continuous(name="Genome position")
  gg
}
