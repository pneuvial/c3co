library(c3co)
library(reshape)
library(tis)
library(ggplot2)
library(matrixStats)
library(mclust)

source("00.functions.R")
source("01.setup.R")
stats <- c("TCN","C1C2", "TCN")
meth <- c("FLLAT","C3CO", "C3CO")
pathFig <- R.utils::Arguments$getWritablePath("Figures")
subClones <- readRDS(file.path(pathSubClones, sprintf("subclones.rds")))
nbClones <- length(subClones)
nbSimu <- length(readRDS(file.path(pathDat, sprintf("simu.rds"))))

###########################################
### Plot latent Profiles
###########################################
message("Sample only 1/10 observations to reduce the size of the figure")
i <- sort(sample(1:len, size = len/10))
df.Sub <- do.call(rbind, lapply(subClones, function(ss) ss[i,]))
df.Sub$feat <- as.factor(rep(1:nbClones, each = length(i)))
dfToPlot <- data.frame(val = c(df.Sub$ct, df.Sub$baft), var = factor(rep(c("c","b"), each = nrow(df.Sub)), levels = c("c","b")), feat = rep(df.Sub$feat, times = 2), pos = rep(df.Sub$pos, times = 2))

p <- ggplot(dfToPlot, aes(pos, val)) + geom_point(cex = 0.4, alpha = 0.2, pch = 19) + facet_grid(var ~ feat, scale = "free") + theme_bw() + ylab("") + scale_x_continuous(name = "Genome position (Mb)",breaks = c(0, 10000, 20000))
p

###########################################
### Evaluation on weights
###########################################
## Loss
message("Compute loss between true and estimated weights")

weightsMat <- readRDS(file.path(pathWeights, sprintf("weightsMat.rds")))
weightsArray <- lossW(nbSimu, meth, stats, weightsMat)
dataArrayWeights <- data.frame(melt(weightsArray))
dataArrayWeights$method <- factor(dataArrayWeights$method, levels = c("FLLAT-TCN","C3CO-TCN","C3CO-C1C2"))
gWeights <- ggplot(dataArrayWeights) + geom_boxplot(aes(x = method, y = value, fill = method)) + ylab("Loss") + xlab("") + theme_bw()
gWeights
ggsave(gWeights,filename = sprintf("%s/weightLoss.pdf",pathFig), width = 7, height = 5)

### RandIndex
message("Compute Rand Index to evaluate capacity to recover weights")
randIndexArray <- randIndW(nbSimu, meth, stats, weightsMat)
dataArrayrandIndex <- data.frame(melt(randIndexArray))
dataArrayrandIndex$method <- factor(dataArrayrandIndex$method, levels = c("FLLAT-TCN", "C3CO-TCN", "C3CO-C1C2"))
gRandIndex <- ggplot(dataArrayrandIndex) + geom_boxplot(aes(x = method, y = value,fill = method)) + ylab("Rand Index") + xlab("") + theme_bw()
gRandIndex
ggsave(gRandIndex,filename=sprintf("%s/RandIndex.pdf",pathFig), width = 7, height = 5)


#######################################################################
### Archetypes ROC curves Full resolution
#######################################################################
message("Compute ROC and AUC on Subclones")
tol <- c(seq(from = 0.0, to = 5, length = 20))
tol <- sort(tol, decreasing = TRUE)


bkpsByClones <- lapply(subClones, function(sss){
  L <- nrow(sss)
  start <- 1:(L - 1)
  end <- 2:L
  which(sss$region[end] != sss$region[start])
})

bkpsByClones[[length(subClones) + 1]] <- len


regionsByClones <- lapply(1:length(subClones), function(iii){
  sss <- subClones[[iii]]
  sss$region[c(bkpsByClones[[iii]], len)]
})   
regionsByClones[[length(subClones) + 1]] <- "(1,1)"

AUCs_arch <- computeAUC(nbSimu, meth, stats, tol, subClones, weightsMat, regionsByClones) 
df.auc.arch <- melt(AUCs_arch)
df.auc.arch$AUC <- df.auc.arch$value
df.auc.arch$method <- factor(df.auc.arch$method, levels=c("FLLAT-TCN","C3CO-TCN","C3CO-C1C2"))
gAUC <- ggplot(df.auc.arch, aes(x=method, y=AUC))+geom_boxplot(aes( fill=method))+theme_bw()+xlab("")
gAUC
ggsave(gAUC, filename=file.path(pathFig, "AUCsZ.pdf"), width=10, height=5)

### PVE
message("Compute the PVEs")
PVEArray <- pveEval(nbSimu, meth, stats)
meanPVE <- apply(PVEArray, 2,colMeans, na.rm=TRUE)
sdPVE <- apply(PVEArray, 2,colSds, na.rm=TRUE)
dataPVE <- cbind(melt(meanPVE), melt(sdPVE))
colnames(dataPVE) <- c("Features", "method", "PVE", "Features","method", "sds")
dataPVE$Features <- dataPVE$Features+1
dataPVE$method <- factor(dataPVE$method, levels=c("FLLAT-TCN","C3CO-TCN","C3CO-C1C2"))
pPVE <-  ggplot(dataPVE)+geom_line(aes(x=Features,y=PVE,colour=method, lty=method))
pPVE <- pPVE+geom_ribbon(aes(ymin=PVE-sds, ymax=PVE+sds,x=Features,fill=method ),alpha=0.3)+ scale_x_continuous(limits=c(2, 10), breaks=seq(from=2, to =10, by=1))+theme_bw()+scale_y_continuous(limits=c(0,1))
pPVE
ggsave(pPVE, filename=file.path(pathFig, "PVEs.pdf"), width=7, height=5)





