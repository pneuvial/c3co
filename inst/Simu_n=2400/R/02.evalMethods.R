library(InCaSCN)
library(reshape)
library(tis)
source("R/00.functions.R")
framework <- "realistic"
forceM <- FALSE
stats <- c("TCN","C1C2", "TCN")
meth <- c("FLLAT","InCaSCN", "InCaSCN")

dataAnnotTP <- loadCnRegionData(dataSet="GSE13372", tumorFraction=1)
dataAnnotN <- loadCnRegionData(dataSet="GSE13372", tumorFraction=0)

pathSubClones <- Arguments$getWritablePath(sprintf("simArchData"))
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
pathWeight <- Arguments$getWritablePath(sprintf("weightData"))
n <- 30

C1N <- dataAnnotN$c*(1-2*abs(dataAnnotN$b-1/2))/2
c1Mean <- by(C1N, dataAnnotN$genotype, mean)[2]
C2N <- dataAnnotN$c*(1+2*abs(dataAnnotN$b-1/2))/2
c2Mean <- by(C2N, dataAnnotN$genotype, mean)[2]
cMean <- mean(dataAnnotN$c)

###########################################
### Evaluation on weights
###########################################
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
pathFig <- Arguments$getWritablePath("Figures/Weights")
ggsave(gWeights,filename=sprintf("%s/weightLoss_n=2400.pdf",pathFig, framework), width=7, height=5)

#######################################################################
### Archetypes ROC curves Full resolution
#######################################################################

tol <- c(seq(from=0.0, to=1, length=25))
tol <- sort(tol, decreasing=TRUE)
rocArrayArchFull <- array(dim=c(100, length(stats), 2,length(tol)),dimnames=list(b=1:100, method=sprintf("%s-%s",meth,stats), ROC=c("tp", "fp")))
rocDataPath <- Arguments$getWritablePath("rocDataInCaSCN")
fileROC <- "rocArray,full,arch.rds"
forceROC <- TRUE

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
                    
pathFig <- Arguments$getWritablePath("Figures/ROC_InCaSCN")
ggsave(gplotROCarchFull, filename=file.path(pathFig, "ROC,archetypes,FullRes_n=2400.pdf"), width=7, height=5)

AUCs_arch <- apply(rocArrayArchFull, 2, function(roc){
  apply(roc,1,ComputeAUC, FPRs)
})

df.auc.arch <- melt(AUCs_arch)
df.auc.arch$AUC <- df.auc.arch$value
ggplot(df.auc.arch)+geom_boxplot(aes(x=method, y=AUC, color=method))


#######################################################################
### Estimated Profiles ROC curves full res
#######################################################################



tol <- seq(from=0, to=1, length=20)
tol <- sort(tol, decreasing=TRUE)


rocDataPath <- Arguments$getWritablePath("rocDataInCaSCN")
fileROC <- "rocArray,full,profiles.rds"
forceROC=TRUE
if(!file.exists(file.path(rocDataPath, fileROC))||forceROC){
  rocArrayProfFull <- array(dim=c(100, length(stats), 2,length(tol)),dimnames=list(b=1:100, method=sprintf("%s-%s",meth,stats), ROC=c("tp", "fp")))
  for(b in 1:100){
    print(b)
    for(ss in 1:length(stats)){
      stat <- stats[ss]
      mm <- meth[ss]      
      filename <- sprintf("weight,n=%s,b=%s.rds", n,b)
      M <- readRDS(file.path(pathWeight, filename))
      newM <- t(apply(M, 1, function(mm) {
        ww <- which(mm<10)
        mm[ww] <- 0
        return(mm)
      }))
      newM <- cbind(newM, 100-rowSums(newM))
      
      composition <- apply(newM, 1, function(mm) {
        ww <- which(mm!=0)
        return(ww)
      })
      
      dataBest <- loadDataBest(mm, stat, framework, b)
      
###### altered loci determination
      alteredLoci <- sapply(composition, function(cc){
        alt <- rowSums(sapply(cc, function(c){
          if(c!=length(bkpsByClones)){
            res <- subClones[[c]]$region!="(1,1)"
          }else{
            res <- rep(FALSE, len) 
          }
          return(res)
        }))        
      })
      alteredLociY <- apply(alteredLoci, 2, function(cc){
        wwA <- which(cc!=0)
        ww <- which(cc==0)
        alt <- logical(0)
        alt[wwA] <- TRUE
        alt[ww] <- FALSE
        return(alt)
      })
      if(stat=="C1C2"){
        Y1hat <- dataBest$Y.hat$Y1
        Y2hat <- dataBest$Y.hat$Y2
        bkp <- dataBest$bkp
        start <- c(1,ceiling(bkp))
        end <- c(floor(bkp), len)
        Y1hatFull <- expand(t(Y1hat),start, end)
        Y2hatFull <- expand(t(Y2hat),start, end)
        YhatFull <- expand(Y1hat+Y2hat,start, end)
        SESP <- SESPC1C2(Y1hatFull,Y2hatFull,alteredLociY,ind=1:ncol(alteredLociY), tol, c1Mean, c2Mean)
              
      }else{
        if(mm=="InCaSCN"){
          YTCNhat <- dataBest$Y.hat$Y1
          bkp <- dataBest$bkp
          start <- c(1,ceiling(bkp))
          end <- c(floor(bkp), len)
          YTCNhatFull <- expand(t(YTCNhat),start, end)
        }else{
          YTCNhatFull <- 2*2^(dataBest$Y.hat)
        }
        SESP <- SESPTCN(YTCNhatFull,alteredLociY,ind=1:ncol(alteredLociY), tol, cMean)
      }     
      rocArrayProfFull[b,ss,"tp",] <- unlist(SESP["tp",])
      rocArrayProfFull[b,ss,"fp",] <- unlist(SESP["fp",])
    }
  }
  saveRDS(rocArrayProfFull, file.path(rocDataPath, fileROC))
}else{
  rocDataPath <- Arguments$getWritablePath("rocDataInCaSCN")
  fileROC <- "rocArray,full,profile.rds"
  rocArrayProfFull <- readRDS(file.path(rocDataPath, fileROC))
}

FPRs <-  sort(unique(as.vector(rocArrayProfFull[,,"fp",])))
sumROCarrayProfFull <- apply(rocArrayProfFull, 2, function(roc) rowMeans(apply(roc,1,ComputeROC, FPRs), na.rm=TRUE))
rownames(sumROCarrayProfFull) <- FPRs
dataROCprofFull <- melt(sumROCarrayProfFull)
colnames(dataROCprofFull) <- c("FPR", "method", "TPR")

gplotROCprofFull <- ggplot(dataROCprofFull)+ geom_line(aes(x=FPR,y=TPR, group=method,colour = method, lty=method))+xlim(c(0,1))+geom_abline(slope=1, intercept=0, colour="red", lty=2)+ylim(c(0,1))+theme_bw()
gplotROCprofFull+geom_point(aes(x=FPR,y=TPR, group=method,colour = method, lty=method))

ggsave(gplotROCprofFull, filename=file.path(pathFig, "ROC,profiles,FullRes_n=2400.pdf"), width=7, height=5)


AUCs <- apply(rocArrayProfFull, 2, function(roc){
  apply(roc,1,ComputeAUC, FPRs)
})

df.auc <- melt(AUCs)
df.auc$AUC <- df.auc$value
ggplot(df.auc)+geom_boxplot(aes(x=method, y=AUC, color=method))


####################################
## Evaluation on Y observed
####################################
rocArrayYFull <- array(dim=c(100, 2, 2,length(tol)),dimnames=list(b=1:100, method=sprintf("%s",stats[1:2]), ROC=c("tp", "fp")))
pathWeight <- Arguments$getWritablePath(sprintf("weightData"))
pathSim <- Arguments$getWritablePath(sprintf("simData"))
for(b in 1:100){
  print(b)
  filename <- sprintf("weight,n=%s,b=%s.rds", n,b)
  M <- readRDS(file.path(pathWeight, filename))
  newM <- t(apply(M, 1, function(mm) {
    ww <- which(mm<10)
    mm[ww] <- 0
    return(mm)
  }))
  newM <- cbind(newM, 100-rowSums(newM))
  
  composition <- apply(newM, 1, function(mm) {
    ww <- which(mm!=0)
    return(ww)
  })
  alteredLoci<- sapply(composition, function(cc){
    alt <- rowSums(sapply(cc, function(c){
      if(c!=length(bkpsByClones)){
        res <- subClones[[c]]$region!="(1,1)"
      }else{
        res <- rep(FALSE, len) 
      }
      return(res)
    }))
    
  })
   alteredLociY <- apply(alteredLoci, 2, function(cc){
        wwA <- which(cc!=0)
        ww <- which(cc==0)
        alt <- logical(0)
        alt[wwA] <- TRUE
        alt[ww] <- FALSE
        return(alt)
      })
  dat <- readRDS(sprintf("%s/dat_B=%s.rds",pathSim,b))
  for(ss in 1:2){
    stat <- stats[ss]    
    resSegmentation <- InCaSCN:::segmentData(dat, stat=stat)
    Y <- resSegmentation$Y
    Y1 <- resSegmentation$Y1
    Y2 <- resSegmentation$Y2
    bkp <- resSegmentation$bkp[[1]]
    start <- c(1,ceiling(bkp))
    end <- c(floor(bkp), len)
    if(stat=="C1C2"){
      start <- c(1,ceiling(bkp))
      end <- c(floor(bkp), len)
      Y1hatFull <- expand(Y1,start, end)
      Y2hatFull <- expand(Y2,start, end)
      YhatFull <- expand(Y1hat+Y2hat,start, end)
      SESP <- SESPC1C2(Y1hatFull,Y2hatFull,alteredLociY,ind=1:ncol(alteredLociY), tol, c1Mean, c2Mean)
      
    }else{
      start <- c(1,ceiling(bkp))
      end <- c(floor(bkp), len)
      YTCNhatFull <- expand(Y,start, end)
      SESP <- SESPTCN(YTCNhatFull,alteredLociY,ind=1:ncol(alteredLociY), tol, cMean)
    }  
    rocArrayYFull[b,ss,"tp",] <- unlist(SESP["tp",])
    rocArrayYFull[b,ss,"fp",] <- unlist(SESP["fp",])   
  }
}


AUCsY <- apply(rocArrayYFull, 2, function(roc){
  apply(roc,1,ComputeAUC, FPRs)
})

df.aucY <- melt(AUCsY)
df.aucY$AUC <- df.aucY$value
ggplot(df.aucY)+geom_boxplot(aes(x=method, y=AUC,fill = factor(stat)))

df.auc.tot <- rbind(cbind(df.auc.arch, stat="Z") , cbind(df.auc, stat="Yhat"))
df.aucY <- melt(colMeans(AUCsY))
df.aucY$Y <- factor(rownames(df.aucY))

gAUC <- ggplot(df.auc.tot, aes(x=stat, y=AUC))+geom_boxplot(aes( fill=method))+theme_bw()+geom_hline(data=df.aucY, aes(yintercept=value, lty=Y),show_guide = TRUE)
gAUC
pathFig <- Arguments$getWritablePath("Figures/AUC")
ggsave(gAUC, filename=file.path(pathFig, "AUCs_n=2400.pdf"), width=10, height=7)



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
    pathArch <- Arguments$getWritablePath(sprintf("archetypeData%s_%s_%s", stat, framework, meth[ss]))
    pathMeth <- Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
    filename <- sprintf("archData_B=%s_%s.rds", b, meth[ss])
    file <- file.path(pathMeth,filename)
    if(file.exists(file)){
      res <- readRDS(file)
      
      
      if(mm=="FLLAT"){
        p.list <- 2:15
        pathfllat <- Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
        listRes <- list()
        listRes$nb.arch <- p.list
        pp=2
        file <- sprintf("%s/featureData,p=%s.rds",pathfllat,pp)
        res <- readRDS(file)
        listRes$PVE <- res$PVE$PVEs
        
      }
      if(mm=="InCaSCN"){
        listRes <- list()
        listRes$PVE <- unlist(sapply(res, function (rr) rr$PVE))
        
      }
      PVEArray[b,ss,1:length(listRes$PVE)] <- listRes$PVE
    }
  }
}

meanPVE <- apply(PVEArray, 2,colMeans, na.rm=TRUE)
sdPVE <- apply(PVEArray, 2,colSds, na.rm=TRUE)

dataPVE <- cbind(melt(meanPVE), melt(sdPVE))[,-c(4,5)]
colnames(dataPVE) <- c("Features", "method", "PVE", "sds")
dataPVE$Features <- dataPVE$Features+1
pPVE <-  ggplot(dataPVE)+geom_line(aes(x=Features,y=PVE,colour=method, lty=method))

pPVE <- pPVE+geom_ribbon(aes(ymin=PVE-sds, ymax=PVE+sds,x=Features,fill=method ),alpha=0.3)+ scale_x_continuous(limits=c(2, 15), breaks=seq(from=2, to =15, by=1))+theme_bw()+scale_y_continuous(limits=c(0.3,1))
pathFig <- Arguments$getWritablePath("Figures/PVE")
ggsave(pPVE, filename=file.path(pathFig, "PVEs_n=2400.pdf"), width=7, height=5)

