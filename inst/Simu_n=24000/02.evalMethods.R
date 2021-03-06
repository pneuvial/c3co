library(c3co)
library(reshape)
library(tis)
library(ggplot2)
library(matrixStats)
message("source 00.functions.R")
source("00.functions.R")
framework <- "realistic"
forceM <- FALSE
Bmax <- 10
B <- 1:Bmax
stats <- c("TCN","C1C2", "TCN")
meth <- c("FLLAT","c3co", "c3co")
pathFig <- R.utils::Arguments$getWritablePath("Figures")
message("Parameters of data sets")
dataAnnotTP <- acnr::loadCnRegionData(dataSet = "GSE13372", tumorFraction = 1)
dataAnnotN <- acnr::loadCnRegionData(dataSet = "GSE13372", tumorFraction = 0)

pathSubClones <- R.utils::Arguments$getWritablePath(sprintf("simArchData"))
subClones <- readRDS(file.path(pathSubClones, list.files(pathSubClones)))
len <- nrow(subClones[[1]])

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
pathWeight <- R.utils::Arguments$getWritablePath(sprintf("weightData"))
n <- 30
message("Compute mean of TCN, minor and major copy numbers in initial data")
C1N <- dataAnnotN$c*(1 - 2*abs(dataAnnotN$b - 1/2))/2
c1Mean <- by(C1N, dataAnnotN$genotype, mean)[2]
C2N <- dataAnnotN$c*(1 + 2*abs(dataAnnotN$b - 1/2))/2
c2Mean <- by(C2N, dataAnnotN$genotype, mean)[2]
cMean <- mean(dataAnnotN$c)
###########################################
### Plot latent Profiles
###########################################
message("Sample only 5000 observations to reduce the size of the figure")
i <- sort(sample(1:len, size = 5000))
df.Sub <- do.call(rbind, lapply(subClones, function(ss) ss[i,]))
df.Sub$feat <- as.factor(rep(1:5, each = length(i)))
dfToPlot <- data.frame(val = c(df.Sub$ct, df.Sub$baft), var = factor(rep(c("c","b"), each = nrow(df.Sub)), levels = c("c","b")), Feat = rep(df.Sub$Feat, times = 2), pos = rep(df.Sub$pos, times = 2))
p <- ggplot(dfToPlot, aes(pos, val)) + geom_point(cex = 0.4, alpha = 0.2, pch = 19) + facet_grid(var ~ feat, scale = "free") + theme_bw() + ylab("") + scale_x_continuous(name = "Genome position (Mb)",breaks = c(0, 10000, 20000))
p

###########################################
### Evaluation on weights
###########################################
## Loss
message("Compute loss between true and estimated weights")
weightsArray <- array(dim = c(Bmax, length(stats)),dimnames = list(b = B, method = sprintf("%s-%s",meth,stats)))
for (b in B) {
  print(b)
  filename <- sprintf("weight,n=%s,b=%s.rds", n,b)
  M <- readRDS(file.path(pathWeight, filename))
  WT <- cbind(M, 100 - rowSums(M))/100
  for (ss in 1:length(stats)) {
    stat = stats[ss]
    mm <- meth[ss]
    dataBest <- loadDataBest(mm, stat, framework, b)
    if (!is.null(dataBest)) {
      W <- dataBest@res@W    
      eps <- 0.1
      corrBestW <- apply(WT,2,function(ww){
        apply(W,2,function(wwh){
          sum(abs(ww - wwh) < eps,na.rm = TRUE)/length(wwh)
        }) 
      })
    ind <- apply(corrBestW, 2, which.max)
    West <- round(W[,ind],2)     
    resLossWeights <- sum((WT - West)^2)/(ncol(West)*ncol(West))
    }else{
      resLossWeights <- NA
    }
    weightsArray[b,ss] <- resLossWeights
  }
}

dataArrayWeights <- data.frame(melt(weightsArray))
dataArrayWeights$method <- factor(dataArrayWeights$method, levels = c("FLLAT-TCN","c3co-TCN","c3co-C1C2"))
gWeights <- ggplot(dataArrayWeights) + geom_boxplot(aes(x = method, y = value, fill = method)) + ylab("Loss") + xlab("") + theme_bw()
gWeights
ggsave(gWeights,filename = sprintf("%s/weightLoss_n=24000.pdf",pathFig, framework), width = 7, height = 5)



### RandIndex
message("Compute Rand Index to evaluate capacity to recover weights")
library(mclust)
randIndexArray <- array(dim = c(Bmax, length(stats)),dimnames = list(b = 1:Bmax, method = sprintf("%s-%s",meth,stats)))
pathWeight <- R.utils::Arguments$getWritablePath(sprintf("weightData"))
pathSim <- R.utils::Arguments$getWritablePath(sprintf("simData"))
for (b in B) {
  print(b)
  filename <- sprintf("weight,n=%s,b=%s.rds", n,b)
  WT <- readRDS(file.path(pathWeight, filename))
  
  WT <- cbind(WT, 100 - rowSums(WT))
  clustWT <-  cutree(hclust(dist(WT),method = "ward.D"), ncol(WT))
  for (ss in 1:length(stats)) {
    stat = stats[ss]
    mm <- meth[ss]
    randIndex <- NA
    dataBest <- loadDataBest(mm, stat, framework, b)
    if (!is.null(dataBest)) {
      W <- dataBest@res@W
      clustW <- cutree(hclust(dist(WT), method = "ward.D"), ncol(W))
      randIndex <- adjustedRandIndex(clustWT, clustW)
    }
    randIndexArray[b,ss] <- randIndex
  }
}
names(randIndexArray) <- c("b", "method")
dataArrayrandIndex <- data.frame(melt(randIndexArray))
dataArrayrandIndex$method <- factor(dataArrayrandIndex$method, levels = c("FLLAT-TCN", "c3co-TCN", "c3co-C1C2"))
gRandIndex <- ggplot(dataArrayrandIndex) + geom_boxplot(aes(x = method, y = value,fill = method)) + ylab("Rand Index") + xlab("") + theme_bw()
gRandIndex
ggsave(gRandIndex,filename=sprintf("%s/RandIndex_n=24000.pdf",pathFig, framework), width = 7, height = 5)


#######################################################################
### Archetypes ROC curves Full resolution
#######################################################################
message("Compute ROC and AUC on Subclones")

tol <- c(seq(from = 0.0, to = 1, length = 20))
tol <- sort(tol, decreasing = TRUE)
rocArrayArchFull <- array(dim = c(Bmax, length(stats), 2,length(tol)),dimnames = list(b = 1:Bmax, method = sprintf("%s-%s",meth,stats), ROC = c("tp", "fp")))
AUCs_arch <- array(dim = c(Bmax, length(stats)),dimnames = list(b = 1:Bmax, method = sprintf("%s-%s",meth,stats)))
rocDataPath <- R.utils::Arguments$getWritablePath("rocDatac3co")
fileROC <- "rocArray,full,arch.rds"
forceROC <- FALSE

alteredLociInClones <- sapply(1:(length(regionsByClones) - 1), function(c){
  if (c != length(regionsByClones)){
    res <- subClones[[c]]$region != "(1,1)"
  }else{
    res <- rep(FALSE, len) 
  }
  return(res)
})

if (!file.exists(file.path(rocDataPath, fileROC)) || forceROC) {
  for (b in B) {
    print(b)
    filename <- sprintf("weight,n=%s,b=%s.rds", n,b)
    M <- readRDS(file.path(pathWeight, filename))
    WT <- cbind(M, 100 - rowSums(M))/100    
    for (ss in 1:length(stats)) {
      stat <- stats[ss]
      mm <- meth[ss]      
      dataBest <- loadDataBest(mm,stat, framework, b)
      if (!is.null(dataBest)) {
        message(sprintf("Compute ROC and AUC for method %s, var %s and data set %s" ,mm, stat, b))
        Z <- dataBest@res@Zt$Z
        W <- dataBest@res@W
        eps <- 0.1
        corrBestW <- apply(WT,2,function(zz){
          apply(W,2,function(zzh){
            sum(abs(zz - zzh) < eps,na.rm = TRUE)/length(zzh)
          }) 
        })
        ind <- apply(corrBestW, 2, which.max)
      
        if (stat == "C1C2") {
          Z1 <- dataBest@res@Zt$Z1
          Z2 <- dataBest@res@Zt$Z2
          bkp <- dataBest@bkp[[1]]
          start <- c(1,ceiling(bkp))
          end <- c(floor(bkp), len)
          Z1hatFull <- expand(Z1,start, end)
          Z2hatFull <- expand(Z2,start, end)
          ZhatFull <- expand(Z,start, end)

          SESP <- SESPC1C2(Z1hatFull,Z2hatFull,alteredLociInClones,ind, tol, 1, 1)
        }else{
          if(mm == "FLLAT") {
            ZhatFull <- t(2*2^Z)          
          }else{
            bkp <- dataBest@bkp
            start <- c(1,ceiling(bkp))
            end <- c(floor(bkp), len)            
            ZhatFull <- expand(Z,start, end)
          }
          SESP <- SESPTCN(ZhatFull,alteredLociInClones, ind,tol, 2)
      
        }
        SESP.mat <- rbind(tp=unlist(SESP["tp",]),fp=unlist(SESP["fp",]))
        AUCs_arch[b,ss] <- ComputeAUC(SESP.mat,unlist(SESP["fp",]))
        rocArrayArchFull[b,ss,"tp",] <- unlist(SESP["tp",])
        rocArrayArchFull[b,ss,"fp",] <- unlist(SESP["fp",])
        }
      }
    }
    saveRDS(rocArrayArchFull, file.path(rocDataPath, fileROC))
  }else{
  rocArrayArchFull <- readRDS(file.path(rocDataPath, fileROC))
}

FPRs <-  sort(unique(as.vector(rocArrayArchFull[,,"fp",])))
pathFig <- R.utils::Arguments$getWritablePath("Figures/ROC_c3co")

AUCs_arch<- apply(rocArrayArchFull, 2, function(roc){
  apply(roc,1,ComputeAUC, FPRs)
})
df.auc.arch <- melt(AUCs_arch)
df.auc.arch$AUC <- df.auc.arch$value
df.auc.arch$method <- factor(df.auc.arch$method, levels=c("FLLAT-TCN","c3co-TCN","c3co-C1C2"))
gAUC <- ggplot(df.auc.arch, aes(x=method, y=AUC))+geom_boxplot(aes( fill=method))+theme_bw()+xlab("")
gAUC
ggsave(gAUC, filename=file.path(pathFig, "AUCs_n=24000,Z.pdf"), width=10, height=5)

### PVE
message("Compute the PVEs")
PVEArray <- array(dim=c(Bmax, length(stats), max(p.list)-1),dimnames=list(b=1:Bmax, method=sprintf("%s-%s",meth,stats)))
for(b in 1:Bmax){
  print(b)
  for(ss in 1:length(stats)){
    stat=stats[ss]
    mm=meth[ss]
    print(sprintf("meth=%s, stat=%s", mm, stat))
    pathMeth <- R.utils::Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
    filename <- sprintf("archData_B=%s_%s.rds", b, meth[ss])
    file <- file.path(pathMeth,filename)
    if(file.exists(file)){
      res <- readRDS(file)
      PVEArray[b,ss,1:length(res)] <- sapply(res, function (rr) rr@PVE)
    }
  }
}

meanPVE <- apply(PVEArray, 2,colMeans, na.rm=TRUE)
sdPVE <- apply(PVEArray, 2,colSds, na.rm=TRUE)

## data to plot PVE
dataPVE <- cbind(melt(meanPVE), melt(sdPVE))
colnames(dataPVE) <- c("Features", "method", "PVE", "Features","method", "sds")
dataPVE$Features <- dataPVE$Features+1
dataPVE$method <- factor(dataPVE$method, levels=c("FLLAT-TCN","c3co-TCN","c3co-C1C2"))
pPVE <-  ggplot(dataPVE)+geom_line(aes(x=Features,y=PVE,colour=method, lty=method))

pPVE <- pPVE+geom_ribbon(aes(ymin=PVE-sds, ymax=PVE+sds,x=Features,fill=method ),alpha=0.3)+ scale_x_continuous(limits=c(2, 15), breaks=seq(from=2, to =15, by=1))+theme_bw()+scale_y_continuous(limits=c(0.15,1))
pPVE
ggsave(pPVE, filename=file.path(pathFig, "PVEs_n=24000.pdf"), width=7, height=5)





