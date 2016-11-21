#' Function to plot PVE
#'
#' @export
#' @param resInCaSCN output from InCaSCN.
#' @param bestNbLatent best number of latent profiles.
#' @param ylim a vector that define min and max of y-axis
#' @return PVE curve
pvePlot <- function(resInCaSCN,bestNbLatent=NULL, ylim=c(0,1)){
  PVEs <- sapply(resInCaSCN, function (rr) rr$PVE)
  nb.arch <- sapply(resInCaSCN, function (rr) rr$param$nb.arch)
  df.PVE <- data.frame(PVE=PVEs, nb.arch=nb.arch)
  gg <- ggplot2::ggplot(df.PVE,aes(x=nb.arch, y=PVE))+ geom_line()+geom_point()+theme_bw()+ylim(ylim)+xlab("Number of latent profiles")
  if(!is.null(bestNbLatent)){
   gg <- gg+ggplot2::geom_vline(xintercept=bestNbLatent, lty=2)
  }
  gg  
}

#' Function to plot W
#'
#' @export
#' @param dataBest A list from output of InCaSCN
#' @param rownamesW A vector that contains identification of patients
#' @param col A vector that contains colors for the heatmap
#' @param margins A vector margins 
#' @param posLegend position of the legend as for \code{plot}
#' @param listPheno A Matrix that contains details on phenotype for each patient. Could be location or time point of tumors for example
#' @param colsPheno Matrix that containts colors for each type of variable in phenotype
#' @return Heatmap of W
Wplot <- function(dataBest, rownamesW=NULL, col= colorRampPalette(RColorBrewer::brewer.pal(9, 'GnBu'))(100),margins=c(5,7),posLegend=NA, listPheno=NA, colsPheno=NA, colLegend, labelLegend){
  W <- dataBest$res$W
  rownames(W) <- rownamesW
  res.clust = hclust(dist(W),method="ward.D")
  if(ncol(listPheno)!=ncol(colsPheno)){
    stop("listPheno and colsPheno must have the same number of columns")
  }
  if(!missing(colsPheno)){
    heatmap.3(W, Rowv=as.dendrogram(res.clust), dendrogram="row",  RowSideColors=t(colsPheno), col=col,scale="none", key=TRUE, cexCol=1.5, cexRow=1.5,margins = c(5,10))
    if(!is.na(posLegend)){
      legend(posLegend,legend=labelLegend, fill=colLegend,border=FALSE, bty="n", y.intersp = 1, cex=1)
    }
  }else{
    heatmap.3(W, Rowv=as.dendrogram(res.clust), dendrogram="row", col=col,scale="none", key=TRUE, cexCol=1.5, cexRow=1.5,margins = margins)
  }  
}

#' Function to plot PVE
#'
#' @export
#' @param minMaxPos Matrix that contains min and max position for each chromosome
#' @param chromsomes A vector that contains the focused chromosomes
#' @param var TCN, Minor or Major 
#' @param dataBest A list from output of InCaSCN
#' @return A data frame to plot Latent profiles with ggplot
createZdf <- function(minMaxPos, dataBest, chromosomes, var="TCN"){
  lengthCHR <- sapply(dataBest$bkp, length)
  chrs <- sapply(1:22, function(cc) rep(cc,times=lengthCHR[cc]))
  start <- c(1,cumsum(lengthCHR)+1)
  
  
  df.CHR <- do.call(rbind, lapply(chromosomes, function(cc){
    bb <- c(minMaxPos[cc,"minPos"], dataBest$bkp[[cc]],minMaxPos[cc,"maxPos"]) 
    if(var=="TCN"){
      zz <- rbind(dataBest$res$Z[start[cc],],dataBest$res$Z[start[cc]:(start[cc+1]-1),],dataBest$res$Z[start[cc+1]-1,])
      
    }else if(var=="Minor"){
      zz <- rbind(dataBest$res$Z1[start[cc],],dataBest$res$Z1[start[cc]:(start[cc+1]-1),],dataBest$res$Z1[start[cc+1]-1,])
      
    }else if(var=="Major"){
      zz <- rbind(dataBest$res$Z2[start[cc],],dataBest$res$Z2[start[cc]:(start[cc+1]-1),],dataBest$res$Z2[start[cc+1]-1,])
  arch <- factor(rep(1:ncol(zz), each=length(bb)))
    }else{
      stop("var must be TCN, Minor or Major")
    }
    arch <- factor(rep(1:ncol(zz), each=length(bb))) 
    zz <- c(zz)
    data.frame(position=bb,CN=zz, arch=arch, chr=factor(sprintf("chromosome %s",cc))) 
  }))
  return(df.CHR)
}

#' Function to plot Latent profiles
#'
#' @export
#' @param df data.frame object output from \code{createZdf}
#' @param ylab Label of y-axis
#' @return plot of Latent profiles
Zplot <- function(df, ylab) {
  gg <- ggplot2::ggplot(df)+ggplot2::geom_step(aes(position, CN, group=arch, col=arch,lty=arch), direction="hv", lwd=1)+ggplot2::facet_wrap(~chr, ncol=2, scale="free")+ggplot2::theme_bw()+xlab("Genome position (Mb)")+ ggplot2::labs(colour = "Latent Profile",lty="Latent Profile")+ggplot2::scale_x_continuous(breaks=seq(from=0,to =150, by=20))+ggplot2::scale_y_continuous(name=ylab)
  gg
}