library(c3co)
library(reshape)
library(tis)
library(ggplot2)
library(matrixStats)
source("R/00.functions.R")
framework <- "realistic"
forceM <- FALSE
stats <- c("TCN","C1C2", "TCN")
meth <- c("FLLAT","c3co", "c3co")
pathFig <- R.utils::Arguments$getWritablePath("Figures")

dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE13372", tumorFraction=1)
dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE13372", tumorFraction=0)

pathSubClones <- R.utils::Arguments$getWritablePath(sprintf("simArchData"))
subClones <- readRDS(file.path(pathSubClones, list.files(pathSubClones)))
len <- nrow(subClones[[1]])
bkpsByClones <- lapply(subClones, function(sss){
  L <- nrow(sss)
  start <- 1:(L-1)
  end <- 2:L
  which(sss$region[end]!=sss$region[start])
})

bkpsByClones[[length(subClones)+1]] <- len
regionsByClones <- lapply(1:length(subClones), function(iii){
  sss <- subClones[[iii]]
  sss$region[c(bkpsByClones[[iii]], len)]
})    
regionsByClones[[length(subClones)+1]] <- "(1,1)"
pathWeight <- R.utils::Arguments$getWritablePath(sprintf("weightData"))
n <- 30

C1N <- dataAnnotN$c*(1-2*abs(dataAnnotN$b-1/2))/2
c1Mean <- by(C1N, dataAnnotN$genotype, mean)[2]
C2N <- dataAnnotN$c*(1+2*abs(dataAnnotN$b-1/2))/2
c2Mean <- by(C2N, dataAnnotN$genotype, mean)[2]
cMean <- mean(dataAnnotN$c)
###########################################
### Plot latent Profiles
###########################################

## Sample only 5000 observations to reduce the size of the figure
i <- sort(sample(1:len, size=5000))
df.Sub <- do.call(rbind, lapply(subClones, function (ss) ss[i,]))
df.Sub$Feat <- as.factor(rep(1:5, each=length(i)))
dfToPlot <- data.frame(val=c(df.Sub$ct, df.Sub$baft),var=factor(rep(c("c","b"), each=nrow(df.Sub)), levels=c("c","b")), Feat=rep(df.Sub$Feat, times=2), pos=rep(df.Sub$pos, times=2))
p <- ggplot(dfToPlot, aes(pos, val)) + geom_point(cex=0.4, alpha=0.2, pch=19)+ facet_grid(var ~ Feat, scale="free")+theme_bw()+ylab("")+scale_x_continuous(name="Genome position (Mb)",breaks = c(0, 10000, 20000))
p

###########################################
### Evaluation on weights
###########################################
## Loss
weightsArray <- array(dim=c(100, length(stats)),dimnames=list(b=1:100, method=sprintf("%s-%s",meth,stats)))
for(b in 1:100){
  print(b)
  filename <- sprintf("weight,n=%s,b=%s.rds", n,b)
  M <- readRDS(file.path(pathWeight, filename))
  WT <- cbind(M, 100-rowSums(M))/100
  for(ss in 1:length(stats)){
    stat=stats[ss]
    mm <- meth[ss]
    dataBest <- loadDataBest(mm, stat, framework, b)
    W <- dataBest$W    
    eps <- 0.1
    corrBestW <- apply(WT,2,function(ww){
      apply(W,2,function(wwh){
        sum(abs(ww-wwh)<eps,na.rm=T)/length(wwh)
      }) 
    })
    ind <- apply(corrBestW, 2, which.max)
    West <- round(W[,ind],2)     
    resLossWeights <- sum((WT-West)^2)/(ncol(West)*ncol(West))
    weightsArray[b,ss] <- resLossWeights
  }
}

dataArrayWeights <- data.frame(melt(weightsArray))
gWeights <- ggplot(dataArrayWeights)+geom_boxplot(aes(x=method, y=value,fill=method))+ylab("Loss")+xlab("")+theme_bw()
gWeights
ggsave(gWeights,filename=sprintf("%s/weightLoss_n=24000.pdf",pathFig, framework), width=7, height=5)



### RandIndex
library(mclust)
randIndexArray <- array(dim=c(100, length(stats)),dimnames=list(b=1:100, method=sprintf("%s-%s",meth,stats)))
pathWeight <- R.utils::Arguments$getWritablePath(sprintf("weightData"))
pathSim <- R.utils::Arguments$getWritablePath(sprintf("simData"))
for(b in 1:100){
  print(b)
  filename <- sprintf("weight,n=%s,b=%s.rds", n,b)
  WT <- readRDS(file.path(pathWeight, filename))
  
  WT <- cbind(WT, 100-rowSums(WT))
  clustWT <-  cutree(hclust(dist(WT),method="ward.D"), ncol(WT))
  for(ss in 1:length(stats)){
    stat=stats[ss]
    mm <- meth[ss]
    dataBest <- loadDataBest(mm, stat, framework, b)
    W <- dataBest$W
    clustW <- cutree(hclust(dist(WT),method="ward.D"), ncol(W))
    randIndex <- adjustedRandIndex(clustWT, clustW)
    randIndexArray[b,ss] <- randIndex
  }
}
randIndexArrayFull <- cbind(randIndexArray)
names(randIndexArrayFull) <- c("b", "method")
dataArrayrandIndex <- data.frame(melt(randIndexArrayFull))
dataArrayrandIndex$method <- factor(dataArrayrandIndex$X2, levels=c("FLLAT-TCN", "c3co-TCN", "c3co-C1C2"))
gRandIndex <- ggplot(dataArrayrandIndex)+geom_boxplot(aes(x=method, y=value,fill=method))+ylab("Rand Index")+xlab("")+theme_bw()
gRandIndex
ggsave(gRandIndex,filename=sprintf("%s/RandIndex_n=24000.pdf",pathFig, framework), width=7, height=5)


#######################################################################
### Archetypes ROC curves Full resolution
#######################################################################

tol <- c(seq(from=0.0, to=1, length=20))
tol <- sort(tol, decreasing=TRUE)
rocArrayArchFull <- array(dim=c(100, length(stats), 2,length(tol)),dimnames=list(b=1:100, method=sprintf("%s-%s",meth,stats), ROC=c("tp", "fp")))
AUCs_arch <- array(dim=c(100, length(stats)),dimnames=list(b=1:100, method=sprintf("%s-%s",meth,stats)))
rocDataPath <- R.utils::Arguments$getWritablePath("rocDatac3co")
fileROC <- "rocArray,full,arch.rds"
forceROC <- FALSE

alteredLociInClones <- sapply(1:(length(regionsByClones)-1), function(c){
  if(c!=length(regionsByClones)){
    res <- subClones[[c]]$region!="(1,1)"
  }else{
    res <- rep(FALSE, len) 
  }
  return(res)
})

if(!file.exists(file.path(rocDataPath, fileROC))||forceROC){
  for(b in 1:100){
    print(b)
    filename <- sprintf("weight,n=%s,b=%s.rds", n,b)
    M <- readRDS(file.path(pathWeight, filename))
    WT <- cbind(M, 100-rowSums(M))/100    
    for(ss in 1:length(stats)){
      stat=stats[ss]
      mm <- meth[ss]      
      dataBest <- loadDataBest(mm,stat, framework, b)
      Z <- dataBest$Z
      W <- dataBest$W
      eps <- 0.1
      corrBestW <- apply(WT,2,function(zz){
        apply(W,2,function(zzh){
          sum(abs(zz-zzh)<eps,na.rm=T)/length(zzh)
        }) 
      })
      ind <- apply(corrBestW, 2, which.max)
      
      if(stat=="C1C2"){
        Z1 <- dataBest$Z1
        Z2 <- dataBest$Z2
        bkp <- dataBest$bkp
        start <- c(1,ceiling(bkp))
        end <- c(floor(bkp), len)
        Z1hatFull <- expand(Z1,start, end)
        Z2hatFull <- expand(Z2,start, end)
        ZhatFull <- expand(Z,start, end)

        SESP <- SESPC1C2(Z1hatFull,Z2hatFull,alteredLociInClones,ind, tol, 1, 1)
      }else{
        if(mm=="FLLAT"){
          ZhatFull <- t(2*2^Z)          
        }else{
          bkp <- dataBest$bkp
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
  saveRDS(rocArrayArchFull, file.path(rocDataPath, fileROC))
}else{
  rocArrayArchFull <- readRDS(file.path(rocDataPath, fileROC))
}

FPRs <-  sort(unique(as.vector(rocArrayArchFull[,,"fp",])))
sumROCarrayArchFull <- apply(rocArrayArchFull, 2, function(roc) rowMeans(apply(roc,1,ComputeROC, FPRs), na.rm=TRUE))
rownames(sumROCarrayArchFull) <- FPRs
dataROCArchFull <- melt(sumROCarrayArchFull)
colnames(dataROCArchFull) <- c("FPR", "method", "TPR")

gplotROCarchFull <- ggplot(dataROCArchFull)+ geom_line(aes(x=FPR,y=TPR, group=method,colour = method, lty=method))+xlim(c(0,1))+geom_abline(slope=1, intercept=0, colour="red", lty=2)+theme_bw()

gplotROCarchFull+ylim(c(0,1))+geom_point(aes(x=FPR,y=TPR, group=method,colour = method))
                    
pathFig <- R.utils::Arguments$getWritablePath("Figures/ROC_c3co")
ggsave(gplotROCarchFull, filename=file.path(pathFig, "ROC,archetypes,FullRes_n=24000.pdf"), width=7, height=5)

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
PVEArray <- array(dim=c(100, length(stats), 15),dimnames=list(b=1:100, method=sprintf("%s-%s",meth,stats)))
for(b in 1:100){
  print(b)
  filename <- sprintf("weight,n=%s,b=%s.rds", n,b)
  M <- readRDS(file.path(pathWeight, filename))
  WT <- cbind(M, 100-rowSums(M))/100
  
  for(ss in 1:length(stats)){
    stat=stats[ss]
    mm=meth[ss]
    print(sprintf("meth=%s, stat=%s", mm, stat))
    pathArch <- R.utils::Arguments$getWritablePath(sprintf("archetypeData%s_%s_%s", stat, framework, meth[ss]))
    pathMeth <- R.utils::Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
    filename <- sprintf("archData_B=%s_%s.rds", b, meth[ss])
    file <- file.path(pathMeth,filename)
    if(file.exists(file)){
      res <- readRDS(file)
      
      
      if(mm=="FLLAT"){
        p.list <- 2:15
        pathfllat <- R.utils::Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
        listRes <- list()
        listRes$nb.arch <- p.list
        pp=2
        file <- sprintf("%s/featureData,p=%s.rds",pathfllat,pp)
        res <- readRDS(file)
        listRes$PVE <- res$PVE$PVEs
        
      }
      if(mm=="c3co"){
        listRes <- list()
        listRes$PVE <- unlist(sapply(res, function (rr) rr$PVE))
        
      }
      PVEArray[b,ss,1:length(listRes$PVE)] <- listRes$PVE
    }
  }
}

meanPVE <- na.omit(apply(PVEArray, 2,colMeans, na.rm=TRUE))
sdPVE <- na.omit(apply(PVEArray, 2,colSds, na.rm=TRUE))

## data to plot PVE
dataPVE <- cbind(melt(meanPVE), melt(sdPVE))[,-c(4,5)]
colnames(dataPVE) <- c("Features", "method", "PVE", "sds")
dataPVE$Features <- dataPVE$Features+1
dataPVE$method <- factor(dataPVE$method, levels=c("FLLAT-TCN","c3co-TCN","c3co-C1C2"))
pPVE <-  ggplot(dataPVE)+geom_line(aes(x=Features,y=PVE,colour=method, lty=method))

pPVE <- pPVE+geom_ribbon(aes(ymin=PVE-sds, ymax=PVE+sds,x=Features,fill=method ),alpha=0.3)+ scale_x_continuous(limits=c(2, 15), breaks=seq(from=2, to =15, by=1))+theme_bw()+scale_y_continuous(limits=c(0.15,1))
pPVE
ggsave(pPVE, filename=file.path(pathFig, "PVEs_n=24000.pdf"), width=7, height=5)





