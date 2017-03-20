##########################################################################
# Run c3co
##########################################################################
## We provide the results of the algorithm and the segmentation so if the user only wants to reproduce figure, run only this file

library("c3co")
library("RColorBrewer")

patient <- "RK29"

### Plot PVE
pvePlot(resC3co@fit, ylim=c(0.70,1))
dataBest <- resC3co@fit[[6]]

### Plot W matrix
figPath <- R.utils::Arguments$getWritablePath("fig-GSE47077-RK29")

pdf(file.path(figPath,sprintf("heatmap,GSE47077,patient=%s.pdf",patient)), width=13, height=8)
Wplot(resC3co, idxBest=6, rownamesW= sprintf("R%s",1:nrow(dataBest@W)))
dev.off()

### Plot Latent profiles
path<- (sprintf("data",patient))
minMaxPos <- readRDS(file.path(path, "minMaxposByCHR.rds"))

lengthCHR <- sapply(resC3co@bkp, length)
chrs <- sapply(1:22, function(cc) rep(cc,times=lengthCHR[cc]))

start <- c(1,cumsum(lengthCHR)+1)

ch <- c(9)
df.CHR <- createZdf(resC3co, minMaxPos, chromosomes=ch, var="TCN", idxBest=6)

df.CHRC1 <- createZdf(resC3co, minMaxPos, chromosomes=ch, var="Minor", idxBest=6)
df.CHRC2 <- createZdf(resC3co, minMaxPos, chromosomes=ch, var="Major", idxBest=6)

df.CHR$position <- df.CHR$position/1e6
df.CHRC2$position <- df.CHRC2$position/1e6
df.CHRC1$position <- df.CHRC1$position/1e6

gArchTCN <- Zplot(df.CHR, ylab="TCN", ylim=c(1,3))
gArchC2 <- Zplot(df.CHRC2, ylab="Major", ylim=c(1,3))
gArchC1 <- Zplot(df.CHRC1, ylab="Minor", ylim=c(0,2))

ggsave(gArchTCN, filename=sprintf("%s/archetypes,GSE47077,patient=%s,chr=%s.pdf", figPath, patient,ch), width=7, height=3.5)
ggsave(gArchC1, filename=sprintf("%s/archetypes,GSE47077,patient=%s,C1,chr=%s.pdf", figPath, patient,ch), width=7, height=3.5)                    
ggsave(gArchC2, filename=sprintf("%s/archetypes,GSE47077,patient=%s,C2,chr=%s.pdf", figPath, patient,ch), width=7, height=3.5)



ch <- c(14)
df.CHR <- createZdf( resC3co,minMaxPos, chromosomes=ch, var="TCN", idxBest=6)

df.CHRC1 <- createZdf( resC3co, minMaxPos, chromosomes=ch, var="Minor", idxBest=6)
df.CHRC2 <- createZdf( resC3co, minMaxPos, chromosomes=ch, var="Major", idxBest=6)

df.CHR$position <- df.CHR$position/1e6
df.CHRC2$position <- df.CHRC2$position/1e6
df.CHRC1$position <- df.CHRC1$position/1e6

gArchTCN <- Zplot(df.CHR, ylab="TCN")
gArchC2 <- Zplot(df.CHRC2, ylab="Major")
gArchC1 <- Zplot(df.CHRC1, ylab="Minor")
ggsave(gArchTCN, filename=sprintf("%s/archetypes,GSE47077,patient=%s,chr=%s.pdf", figPath, patient,ch), width=7, height=3.5)
ggsave(gArchC1, filename=sprintf("%s/archetypes,GSE47077,patient=%s,C1,chr=%s.pdf", figPath, patient,ch), width=7, height=3.5)                    
ggsave(gArchC2, filename=sprintf("%s/archetypes,GSE47077,patient=%s,C2,chr=%s.pdf", figPath, patient,ch), width=7, height=3.5)

