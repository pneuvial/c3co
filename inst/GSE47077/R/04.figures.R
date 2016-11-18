##########################################################################
# Run InCaSCN
##########################################################################
## We provide the results of the algorithm and the segmentation so if the user only wants to reproduce figure, run only this file

library(InCaSCN)
library(RColorBrewer)

patient <- "RK29"
path <- "data"
#dat <- readRDS(sprintf("%s/dat-%s.rds", path ,patient))
output.dir <- Arguments$getWritablePath(sprintf("resultsInCaSCN-%s",patient))
resInC <- list.files(output.dir, pattern="featureData")
resInCaSCN <- lapply(file.path(output.dir, resInC), readRDS)

### Plot PVE
pvePlot(resInCaSCN, ylim=c(0.70,1))
dataBest <- resInCaSCN[[6]]

### Plot W matrix
pathFig <- Arguments$getWritablePath("fig-GSE47077-RK29")

pdf(file.path(pathFig,sprintf("heatmap,GSE47077,patient=%s.pdf",patient)), width=13, height=8)
Wplot(dataBest, rownamesW= gsub(".*-","", names(dat)))
dev.off()

### Plot Latent profiles
minMaxPos <- readRDS(file.path(path, "minMaxposByCHR.rds")

lengthCHR <- sapply(dataBest$bkp, length)
chrs <- sapply(1:22, function(cc) rep(cc,times=lengthCHR[cc]))

start <- c(1,cumsum(lengthCHR)+1)

ch <- c(9)
df.CHR <- createZdf(minMaxPos, dataBest, start, chromosomes=ch)

df.CHRC1 <- createZdf(minMaxPos, dataBest, start, chromosomes=ch, var="Minor")
df.CHRC2 <- createZdf(minMaxPos, dataBest, start, chromosomes=ch, var="Major")

df.CHR$position <- df.CHR$position/1e6
df.CHRC2$position <- df.CHRC2$position/1e6
df.CHRC1$position <- df.CHRC1$position/1e6

gArchTCN <- Zplot(df.CHR, ylab="TCN")
gArchC2 <- Zplot(df.CHRC2, ylab="Major")
gArchC1 <- Zplot(df.CHRC1, ylab="Minor")

ggsave(gArchTCN, filename=sprintf("%s/archetypes,GSE47077,patient=%s,chr=%s.pdf", pathFig, patient,ch), width=7, height=3.5)
ggsave(gArchC1, filename=sprintf("%s/archetypes,GSE47077,patient=%s,C1,chr=%s.pdf", pathFig, patient,ch), width=7, height=3.5)                    
ggsave(gArchC2, filename=sprintf("%s/archetypes,GSE47077,patient=%s,C2,chr=%s.pdf", pathFig, patient,ch), width=7, height=3.5)



ch <- c(14)
df.CHR <- createZdf(minMaxPos, dataBest, start, chromosomes=ch)

df.CHRC1 <- createZdf(minMaxPos, dataBest, start, chromosomes=ch, var="Minor")
df.CHRC2 <- createZdf(minMaxPos, dataBest, start, chromosomes=ch, var="Major")

df.CHR$position <- df.CHR$position/1e6
df.CHRC2$position <- df.CHRC2$position/1e6
df.CHRC1$position <- df.CHRC1$position/1e6

gArchTCN <- Zplot(df.CHR, ylab="TCN")
gArchC2 <- Zplot(df.CHRC2, ylab="Major")
gArchC1 <- Zplot(df.CHRC1, ylab="Minor")
ggsave(gArchTCN, filename=sprintf("%s/archetypes,GSE47077,patient=%s,chr=%s.pdf", pathFig, patient,ch), width=7, height=3.5)
ggsave(gArchC1, filename=sprintf("%s/archetypes,GSE47077,patient=%s,C1,chr=%s.pdf", pathFig, patient,ch), width=7, height=3.5)                    
ggsave(gArchC2, filename=sprintf("%s/archetypes,GSE47077,patient=%s,C2,chr=%s.pdf", pathFig, patient,ch), width=7, height=3.5)

