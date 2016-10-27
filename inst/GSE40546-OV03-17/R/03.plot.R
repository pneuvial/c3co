library(ggplot2)
library("GMD")
library(GEOquery)
library(RColorBrewer)
library(aroma.cn)
source("R/heatmap.3.R")
source("R/01.loadData.R")
minMaxPos <- do.call(rbind, lapply(1:22, function(chromosome){
  chrTag <- sprintf("Chr%s", chromosome)
  units <- getUnitsOnChromosome(ugp, chromosome=chromosome)
  platform <- getPlatform(ugp)
  unf <- getUnitNamesFile(ugp)
  unitNames <- getUnitNames(unf, units=units)
  ## Identify SNP units
  snpPattern <- "^SNP"
  keep <- (regexpr(snpPattern, unitNames) != -1)
  units <- units[keep]
  pos <- getPositions(ugp, units=units)
  return(list(chr=chromosome, minPos=min(pos), maxPos=max(pos)))
}))

minMaxPos

resInC <- list.files(output.dir, pattern="featureData")
resInCaSCN <- lapply(file.path(output.dir, resInC), readRDS)

PVEs <- sapply(resInCaSCN, function (rr) rr$PVE)
nb.arch <- sapply(resInCaSCN, function (rr) rr$param$nb.arch)
df.PVE <- data.frame(PVE=PVEs, nb.arch=nb.arch)
gg <- ggplot(df.PVE,aes(x=nb.arch, y=PVE))+ geom_line()+geom_point()+theme_bw()+ylim(c(0.9,1))+xlab("Number of latent profiles")+geom_vline(xintercept=5, lty=2)
gg
figPath <- Arguments$getWritablePath("Figures-v2")
ggsave(gg, filename=sprintf("%s/PVE,Markowetz.pdf", figPath), width=5, height=3)


idxBestList <- 4
dataBest <- resInCaSCN[[idxBestList]]
res.clust = hclust(dist(dataBest$res$W),method="ward.D")
subInfoClin <- identity[gsm,]
loc <- factor(subInfoClin$characteristics_ch1.3, levels=unique(subInfoClin$characteristics_ch1.3))
timePoint <- factor(subInfoClin$characteristics_ch1.4, levels=unique(subInfoClin$characteristics_ch1.4))
col = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
col1 <- brewer.pal(nlevels(loc), "Set1")
col2 <- brewer.pal(nlevels(timePoint), "Spectral")
W <- round(dataBest$res$W,2)
rownames(W) <- gsub("sample id: ","",subInfoClin$characteristics_ch1.1)
figPath <- Arguments$getWritablePath("figuresGSE40546,C1C2")

pdf(sprintf("%s/heatmap,GSE40546,patient=%s.pdf", figPath, patient), width=13, height=8)
colSide <- t(cbind(loc=col1[loc],timePoint=col2[timePoint]))
rownames(colSide) <- c("Tissue", "Time point")
heatmap.3(W, Rowv=as.dendrogram(res.clust), dendrogram="row",  RowSideColors=colSide, col=col,scale="none", key=TRUE, cexCol=1.5, cexRow=1.5,margins = c(5,10))
legend("topright",legend=c(levels(loc), levels(timePoint)), fill=c(col1[-(nlevels(loc)+1)],col2[-(nlevels(timePoint)+1)]),border=FALSE, bty="n", y.intersp = 1, cex=1)
dev.off()
  
lengthCHR <- sapply(dataBest$bkp, length)
chrs <- rep(1:22, times=lengthCHR)

start <- c(1,cumsum(lengthCHR)+1)
df.CHR <- do.call(rbind, lapply(seq(along=lengthCHR), function(cc){
  bb <- c(minMaxPos[cc,"minPos"][[1]], dataBest$bkp[[cc]],minMaxPos[cc,"maxPos"][[1]]) 
  zz <- rbind(dataBest$res$Z[start[cc],],dataBest$res$Z[start[cc]:(start[cc+1]-1),],dataBest$res$Z[start[cc+1]-1,])
  arch <- factor(rep(1:ncol(zz), each=length(bb)))
  
  zz <- c(zz)
  data.frame(bkp=bb,Z=zz, arch=arch, chr=factor(sprintf("chromosome %s",cc))) 
}))

gArch <- ggplot(df.CHR)+geom_step(aes(bkp, Z, group=arch, col=arch,lty=arch), direction="hv", lwd=2)+facet_wrap(~chr, ncol=2, scale="free")+theme_bw()
gArch

ggsave(gArch, filename=sprintf("%s/archetypes,GSE40546,patient=%s.pdf", figPath, patient), width=15, height=20)



df.CHRC1 <- do.call(rbind, lapply(seq(along=lengthCHR), function(cc){
  bb <- c(minMaxPos[cc,"minPos"][[1]], dataBest$bkp[[cc]],minMaxPos[cc,"maxPos"][[1]]) 
  zz <- rbind(dataBest$res$Z1[start[cc],],dataBest$res$Z1[start[cc]:(start[cc+1]-1),],dataBest$res$Z1[start[cc+1]-1,])
  arch <- factor(rep(1:ncol(zz), each=length(bb)))
  
  zz <- c(zz)
  data.frame(bkp=bb,Z=zz, arch=arch, chr=factor(sprintf("chromosome %s",cc))) 
}))

gArch <- ggplot(df.CHRC1)+geom_step(aes(bkp, Z, group=arch, col=arch,lty=arch), direction="hv", lwd=2)+facet_wrap(~chr, ncol=2, scale="free")+theme_bw()
gArch

ggsave(gArch, filename=sprintf("%s/archetypes,GSE40546,patient=%s,C1.pdf", figPath, patient), width=15, height=20)


df.CHRC2 <- do.call(rbind, lapply(seq(along=lengthCHR), function(cc){
  bb <- c(minMaxPos[cc,"minPos"][[1]], dataBest$bkp[[cc]],minMaxPos[cc,"maxPos"][[1]]) 
  zz <- rbind(dataBest$res$Z2[start[cc],],dataBest$res$Z2[start[cc]:(start[cc+1]-1),],dataBest$res$Z2[start[cc+1]-1,])
  arch <- factor(rep(1:ncol(zz), each=length(bb)))
  
  zz <- c(zz)
  data.frame(bkp=bb,Z=zz, arch=arch, chr=factor(sprintf("chromosome %s",cc))) 
}))

gArch <- ggplot(df.CHRC2)+geom_step(aes(bkp, Z, group=arch, col=arch,lty=arch), direction="hv", lwd=2)+facet_wrap(~chr, ncol=2, scale="free")+theme_bw()
gArch


ggsave(gArch, filename=sprintf("%s/archetypes,GSE40546,patient=%s,C2.pdf", figPath, patient), width=15, height=20)

df.CHRC1C2 <- cbind(df.CHRC2,df.CHRC1)
colnames(df.CHRC1C2) <- c("bkp", "C2","arch", "chr","bkp", "C1","arch", "chr")

gC1C2 <- ggplot(df.CHRC1C2)+geom_line(aes(C1, C2, group=arch, col=arch), direction="hv", lwd=1)+facet_wrap(~chr, ncol=3, scale="free_y")+theme_bw()
gC1C2
ggsave(gC1C2, filename=sprintf("%s/archetypes,GSE40546,patient=%s,C1C2.pdf", figPath, patient), width=15, height=20)




ch <- 9
df.CHR <- do.call(rbind, lapply(ch, function(cc){
  bb <- c(minMaxPos[cc,"minPos"][[1]], dataBest$bkp[[cc]],minMaxPos[cc,"maxPos"][[1]]) 
  zz <- rbind(dataBest$res$Z[start[cc],],dataBest$res$Z[start[cc]:(start[cc+1]-1),],dataBest$res$Z[start[cc+1]-1,])
  arch <- factor(rep(1:ncol(zz), each=length(bb))) 
  zz <- c(zz)
  data.frame(position=bb,TCN=zz, arch=arch, chr=factor(sprintf("chromosome %s",cc))) 
}))

df.CHRC1 <- do.call(rbind, lapply(ch, function(cc){
  bb <- c(minMaxPos[cc,"minPos"][[1]], dataBest$bkp[[cc]],minMaxPos[cc,"maxPos"][[1]]) 
  zz <- rbind(dataBest$res$Z1[start[cc],],dataBest$res$Z1[start[cc]:(start[cc+1]-1),],dataBest$res$Z1[start[cc+1]-1,])
  arch <- factor(rep(1:ncol(zz), each=length(bb))) 
  zz <- c(zz)
  data.frame(position=bb,C1=zz, arch=arch, chr=factor(sprintf("chromosome %s",cc))) 
}))
df.CHRC2 <- do.call(rbind, lapply(ch, function(cc){
  bb <- c(minMaxPos[cc,"minPos"][[1]], dataBest$bkp[[cc]],minMaxPos[cc,"maxPos"][[1]]) 
  zz <- rbind(dataBest$res$Z2[start[cc],],dataBest$res$Z2[start[cc]:(start[cc+1]-1),],dataBest$res$Z2[start[cc+1]-1,])
  arch <- factor(rep(1:ncol(zz), each=length(bb)))
  zz <- c(zz)
  data.frame(position=bb,C2=zz, arch=arch, chr=factor(sprintf("chromosome %s",cc))) 
}))

library(gridExtra)
df.CHR$position <- df.CHR$position/1e6
df.CHRC2$position <- df.CHRC2$position/1e6
df.CHRC1$position <- df.CHRC1$position/1e6

gArchTCN <- ggplot(df.CHR)+geom_step(aes(position, TCN, group=arch, col=arch,lty=arch), direction="hv", lwd=2)+facet_wrap(~chr, ncol=2, scale="free")+theme_bw()+xlab("Genome position (Mb)")+ labs(colour = "Latent Profile",lty="Latent Profile")+scale_x_continuous(breaks=seq(from=0,to =150, by=20))
gArchC2 <- ggplot(df.CHRC2)+geom_step(aes(position, C2, group=arch, col=arch,lty=arch), direction="hv", lwd=2)+facet_wrap(~chr, ncol=2, scale="free")+theme_bw()+xlab("Genome position (Mb)")+ylab("C2")+ labs(colour = "Latent Profile",lty="Latent Profile")+scale_x_continuous(breaks=seq(from=0,to =150, by=20))
gArchC1 <- ggplot(df.CHRC1)+geom_step(aes(position, C1, group=arch, col=arch,lty=arch), direction="hv", lwd=2)+facet_wrap(~chr, ncol=2, scale="free")+theme_bw()+xlab("Genome position (Mb)")+ylab("C1")+ labs(colour = "Latent Profile",lty="Latent Profile")+scale_x_continuous(breaks=seq(from=0,to =150, by=20))
                                                                                                                                                                                                                                                         
grid.arrange(gArchTCN,gArchC1,gArchC2)



#### Save
ggsave(gArchTCN, filename=sprintf("%s/archetypes,GSE40546,patient=%s,chr=%s.pdf", figPath, patient,ch), width=7, height=3.5)
ggsave(gArchC1, filename=sprintf("%s/archetypes,GSE40546,patient=%s,C1,chr=%s.pdf", figPath, patient,ch), width=7, height=3.5)                    
ggsave(gArchC2, filename=sprintf("%s/archetypes,GSE40546,patient=%s,C2,chr=%s.pdf", figPath, patient,ch), width=7, height=3.5)
