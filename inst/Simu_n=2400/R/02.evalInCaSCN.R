library(InCaSCN)
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


alteredLoci <- sapply(1:length(regionsByClones), function(c){
  if(c!=length(regionsByClones)){
    res <- subClones[[c]]$region!="(1,1)"
  }else{
    res <- rep(FALSE, len) 
  }
  return(res)
})

C1N <- dataAnnotN$c*(1-2*abs(dataAnnotN$b-1/2))/2
c1Mean <- by(C1N, dataAnnotN$genotype, mean)[2]
C2N <- dataAnnotN$c*(1+2*abs(dataAnnotN$b-1/2))/2
c2Mean <- by(C2N, dataAnnotN$genotype, mean)[2]
cMean <- mean(dataAnnotN$c)


weightsArray <- array(dim=c(100, length(stats)),dimnames=list(b=1:100, method=sprintf("%s-%s",meth,stats)))
for(b in 1:100){
  print(b)
  filename <- sprintf("weight,n=%s,b=%s.rds", n,b)
  M <- readRDS(file.path(pathWeight, filename))
  WT <- cbind(M, 100-rowSums(M))/100
  for(ss in 1:length(stats)){
    stat=stats[ss]
    pathArch <- Arguments$getWritablePath(sprintf("archetypeData%s_%s_%s", stat, framework, meth[ss]))
    pathMeth <- Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
    filename <- sprintf("archData_B=%s_%s.rds", b, meth[ss])
    file <- file.path(pathMeth,filename)
    if(file.exists(file)){
      res <- readRDS(file)
      
      if(meth[ss]=="FLLAT"){
        p.list <- 2:15
        pathfllat <- Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
        listRes <- list()
        listRes$nb.arch <- p.list
        pp=2
        file <- sprintf("%s/featureData,p=%s.rds",pathfllat,pp)
        res <- readRDS(file)
        listRes$PVE <- res$PVE$PVEs
        pBest <- min(which(diff(listRes$PVE)<1e-3))
        listRes$res <- lapply(p.list, function(pp){
          file <- sprintf("%s/featureData,p=%s.rds",pathfllat,pp)
          res <- readRDS(file)$res
        })
        }
      if(meth[ss]=="InCaSCN"){
        listRes <- list()
        listRes$PVE <- unlist(sapply(res, function (rr) rr$PVE))
        listRes$nb.arch <- unlist(sapply(res, function (rr) rr$param$nb.arch))
        listRes$res <- lapply(res, function (rr) rr$res)
        pBest <- min(c(which(diff(listRes$PVE)<1e-3),length(listRes$PVE) ))
      }
      dataBest <- listRes$res[[pBest]]
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
}

  
library(reshape)
dataArrayWeights <- data.frame(melt(weightsArray))
gWeights <- ggplot(dataArrayWeights)+geom_boxplot(aes(x=method, y=value,fill=method))+ylab("Loss")+xlab("")+theme_bw()
gWeights
pathFig <- Arguments$getWritablePath("Figures/Weights")
ggsave(gWeights,filename=sprintf("%s/weightLoss.pdf",pathFig, framework), width=7, height=5)#######################################################################
### ARchetypes ROC curves
#######################################################################
ComputeROC <- function(roc, FPRs) {
  TPRs <- sapply(FPRs, FUN=function(fp) {
      ww <- which(roc["fp",]<=fp)
      max(roc["tp",ww])
  })
  TPRs[is.infinite(TPRs)] <- 0
   return(TPRs)
}


tol <- c(seq(from=0.0, to=2, length=25))
tol <- sort(tol, decreasing=TRUE)

#######################################################################
### Archetypes ROC curves Full resolution
#######################################################################
rocArrayArchFull <- array(dim=c(100, length(stats), 2,length(tol)),dimnames=list(b=1:100, method=sprintf("%s-%s",meth,stats), ROC=c("tp", "fp")))
rocDataPath <- Arguments$getWritablePath("rocDataInCaSCN")
fileROC <- "rocArray,full,arch.rds"
forceROC=T

if(!file.exists(file.path(rocDataPath, fileROC))||forceROC){
  for(b in 1:100){
    
    print(b)
    filename <- sprintf("weight,n=%s,b=%s.rds", n,b)
    M <- readRDS(file.path(pathWeight, filename))
    WT <- cbind(M, 100-rowSums(M))/100
    
    for(ss in 1:length(stats)){
      stat=stats[ss]
      mm <- meth[ss]
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
          pBest <- min(which(diff(listRes$PVE)<1e-3))
          listRes$res <- lapply(p.list, function(pp){
            file <- sprintf("%s/featureData,p=%s.rds",pathfllat,pp)
            res <- readRDS(file)$res
          })
        }
        if(mm=="InCaSCN"){
          listRes <- list()
          listRes$PVE <- unlist(sapply(res, function (rr) rr$PVE))
          listRes$nb.arch <- unlist(sapply(res, function (rr) rr$param$nb.arch))
          listRes$res <- lapply(res, function (rr) rr$res)
        }
        pBest <- min(c(which(diff(listRes$PVE)<1e-3),length(listRes$PVE) ))
        dataBest <- listRes$res[[pBest]]
        Z <- dataBest$Z
        W <- dataBest$W  
        if(stat=="C1C2"){
          Z1 <- dataBest$Z1
          Z2 <- dataBest$Z2
        }
        eps <- 0.1
        if(stat=="C1C2"){
          bkp <- res[[pBest]]$bkp
          start <- c(1,bkp[[1]]+1)
          end <- c(bkp[[1]], len)
          
          Z1hatFull <- t(apply(Z1, 2, function(yy){
            y1 <- unlist(sapply(seq(along=start), function(ii){
              rep(yy[ii], times=(end[ii]-start[ii]+1))
            }))
          }))
          Z2hatFull <- t(apply(Z2, 2, function(yy){
            y1 <- unlist(sapply(seq(along=start), function(ii){
              rep(yy[ii], times=(end[ii]-start[ii]+1))
            }))
          }))
          estTCN <- t(apply(Z, 2, function(yy){
            y1 <- unlist(sapply(seq(along=start), function(ii){
              rep(yy[ii], times=(end[ii]-start[ii]+1))
            }))
          }))
        }else{
          if(meth[[ss]]=="FLLAT"){
            trueTCN <- sapply(subClones, function(ss){
              c <- ss$ct
            })
            estTCN <- t(2*2^Z)
           
          }else{
            bkp <- res[[pBest]]$bkp
            start <- c(1,ceiling(bkp[[1]]))
            end <- c(floor(bkp[[1]]), len)            
            estTCN <- t(apply(Z, 2, function(yy){
              y1 <- unlist(sapply(seq(along=start), function(ii){
                rep(yy[ii], times=(end[ii]-start[ii]+1))
              }))
            }))       
          }
        }      
        corrBestW <- apply(WT,2,function(zz){
          apply(W,2,function(zzh){
            sum(abs(zz-zzh)<eps,na.rm=T)/length(zzh)
          }) 
        })
        ind <- apply(corrBestW, 2, which.max)
        SESP <- sapply(tol, function(tt){
          getTPTN <- rowSums(sapply(1:ncol(alteredLoci), function(j){
            k <- ind[j]
            regJ <- alteredLoci[,j]
            if(stat=="C1C2"){
              zz1 <- Z1hatFull[k,]
              zz2 <- Z2hatFull[k,]
              wwA <- which(regJ)                
              TP <- sum((abs(zz1[wwA]-c1Mean)>tt | abs(zz2[wwA]-c2Mean)>tt), na.rm=TRUE)
              ww <- which(!regJ)      
              FP <- sum((abs(zz1[ww]-c1Mean)>tt | abs(zz2[ww]-c2Mean)>tt))
            }else{
              zz1 <- estTCN[k,]
              wwA <- which(regJ)                
              TP <- sum(abs(zz1[wwA]-2)>tt, na.rm=TRUE)
              ww <- which(!regJ)      
              FP <- sum(abs(zz1[ww]-2)>tt, na.rm=TRUE)
            }       
            return(c(FP=FP,TP=TP, pos=length(wwA), neg=length(ww)))
          }), na.rm=T)
          se <- getTPTN["TP"]/getTPTN["pos"]
          sp <- getTPTN["FP"]/getTPTN["neg"]
          return(list(tp=se,fp=sp))
        })
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
sumROCarrayArchFull <- apply(rocArrayArchFull, 2, function(roc) rowMeans(apply(roc,1,ComputeROC, FPRs), na.rm=TRUE))
rownames(sumROCarrayArchFull) <- FPRs
dataROCArchFull <- melt(sumROCarrayArchFull)
colnames(dataROCArchFull) <- c("FPR", "method", "TPR")

gplotROCarchFull <- ggplot(dataROCArchFull)+ geom_line(aes(x=FPR,y=TPR, group=method,colour = method, lty=method))+xlim(c(0,1))+geom_abline(slope=1, intercept=0, colour="red", lty=2)+theme_bw()

gplotROCarchFull+ylim(c(0,1))+geom_point(aes(x=FPR,y=TPR, group=method,colour = method))
                    
pathFig <- Arguments$getWritablePath("Figures/ROC_InCaSCN")
ggsave(gplotROCarchFull, filename=file.path(pathFig, "ROC,archetypes,FullRes.pdf"), width=7, height=5)

library(tis)
library(reshape)
ComputeAUC <- function(roc, FPRs) {
    
  TPRs <- sapply(FPRs, FUN=function(fp) {
    ww <- which(roc["fp",]<=fp)
    max(roc["tp",ww])
  })
  TPRs[is.infinite(TPRs)] <- NaN
  auc <- sum(lintegrate(FPRs, TPRs, xint=FPRs)) 
  return(auc)
}

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
forceROC=F
if(!file.exists(file.path(rocDataPath, fileROC))||forceROC){
  rocArrayProfFull <- array(dim=c(100, length(stats), 2,length(tol)),dimnames=list(b=1:100, method=sprintf("%s-%s",meth,stats), ROC=c("tp", "fp")))
  for(b in 1:100){
    print(b)
    for(ss in 1:length(stats)){
      stat=stats[ss]
      mm <- meth[ss]
      print(sprintf("meth=%s, stat=%s", mm, stat))
      pathArch <- Arguments$getWritablePath(sprintf("archetypeData%s_%s_%s", stat, framework, meth[ss]))
      pathMeth <- Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
      filename <- sprintf("archData_B=%s_%s.rds", b, meth[ss])
      file <- file.path(pathMeth,filename)
      
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
          pBest <- min(which(diff(listRes$PVE)<1e-3))
          listRes$res <- lapply(p.list, function(pp){
            file <- sprintf("%s/featureData,p=%s.rds",pathfllat,pp)
            res <- readRDS(file)$res
          })
        }
        if(mm=="InCaSCN"){
          listRes <- list()
          listRes$PVE <- unlist(sapply(res, function (rr) rr$PVE))
          listRes$nb.arch <- unlist(sapply(res, function (rr) rr$param$nb.arch))
          listRes$res <- lapply(res, function (rr) rr$res)
        }
         pBest <- min(c(which(diff(listRes$PVE)<1e-3),length(listRes$PVE) ))
        dataBest <- listRes$res[[pBest]]
          
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
        
        if(stat=="C1C2"){
          Y1hat <- dataBest$Y.hat$Y1
          Y2hat <- dataBest$Y.hat$Y2
          bkp <- res[[pBest]]$bkp
          start <- c(1,ceiling(bkp[[1]]))
          end <- c(floor(bkp[[1]]), len)
          
          Y1hatFull <- t(apply(Y1hat, 1, function(yy){
            y1 <- unlist(sapply(seq(along=start), function(ii){
              rep(yy[ii], times=(end[ii]-start[ii]+1))
            }))
          }))
          Y2hatFull <- t(apply(Y2hat, 1, function(yy){
            y2 <- unlist(sapply(seq(along=start), function(ii){
              rep(yy[ii], times=(end[ii]-start[ii]+1))
            }))
          }))
          
        }else{
          if(mm=="InCaSCN"){
            YTCNhat <- dataBest$Y.hat$Y1
            bkp <- res[[pBest]]$bkp
            start <- c(1,ceiling(bkp[[1]]))
            end <- c(floor(bkp[[1]]), len)
            YTCNhatFull <- t(apply(YTCNhat, 1, function(yy){
              yTCN <- unlist(sapply(seq(along=start), function(ii){
                rep(yy[ii], times=(end[ii]-start[ii]+1))
              }))
            }))
          }else{
            YTCNhatFull <- 2*2^(dataBest$Y.hat)
          }
        }     
        SESP <- sapply(tol, function(tt){
          getTPTN <- rowSums(sapply(1:30, function(j){
            if(stat=="C1C2"){
              wwA <- which(alteredLoci[,j]!=0)
              pos <- (abs(Y1hatFull[j,]-c1Mean)>tt | abs(Y2hatFull[j,]-c2Mean)>tt)
              TP <- sum(pos[wwA], na.rm=TRUE)
              ww <- which(alteredLoci[,j]==0)     
              FP <- sum(pos[ww])
            }else{
              zz1 <- YTCNhatFull[j,]
              wwA <- which(alteredLoci[,j]!=0)
              pos <- abs(zz1-(cMean))>tt
              TP <- sum(pos[wwA], na.rm=TRUE)
              ww <- which(alteredLoci[,j]==0)     
              FP <- sum(pos[ww], na.rm=TRUE)
            }       
            return(c(FP=FP,TP=TP, pos=length(wwA), neg=length(ww)))
          }), na.rm=T)
          se <- getTPTN["TP"]/getTPTN["pos"]
          sp <- getTPTN["FP"]/getTPTN["neg"]
          return(list(tp=se,fp=sp))
        })
        rocArrayProfFull[b,ss,"tp",] <- unlist(SESP["tp",])
        rocArrayProfFull[b,ss,"fp",] <- unlist(SESP["fp",])
      }
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

ggsave(gplotROCprofFull, filename=file.path(pathFig, "ROC,profiles,FullRes.pdf"), width=7, height=5)


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
  dat <- readRDS(sprintf("%s/dat_B=%s.rds",pathSim,b))
  dat <- lapply(dat, function (dd) {
    dd$chr <- 1
    dd$pos <- 1:nrow(dd)
    return(dd)
  })
  for(ss in 1:2){
    stat <- stats[ss]    
    resSegmentation <- segmentData(dat, stat=stat)
    Y <- resSegmentation$Y
    Y1 <- resSegmentation$Y1
    Y2 <- resSegmentation$Y2
    bkp <- resSegmentation$bkp[[1]]
    start <- c(1,ceiling(bkp))
    end <- c(floor(bkp), len)

    Y1hatFull <- t(apply(t(Y1), 1, function(yy){
            yTCN <- unlist(sapply(seq(along=start), function(ii){
              rep(yy[ii], times=(end[ii]-start[ii]+1))
            }))
          }))
    Y2hatFull <- t(apply(t(Y2), 1, function(yy){
      yTCN <- unlist(sapply(seq(along=start), function(ii){
        rep(yy[ii], times=(end[ii]-start[ii]+1))
      }))
    }))
    
    YTCNhatFull <- t(apply(t(Y), 1, function(yy){
      yTCN <- unlist(sapply(seq(along=start), function(ii){
        rep(yy[ii], times=(end[ii]-start[ii]+1))
      }))
    }))

    
    
    SESP <- sapply(tol, function(tt){
      getTPTN <- rowSums(sapply(1:30, function(j){
        regJ <- alteredLoci[,j]
        if(stat=="C1C2"){         
          wwA <- which(regJ!=0)
          pos <- (abs(Y1hatFull[j,]-c1Mean)>tt | abs(Y2hatFull[j,]-c2Mean)>tt)
          TP <- sum(pos[wwA], na.rm=TRUE)
          ww <- which(regJ==0)     
          FP <- sum(pos[ww])
        }else{
          
          zz1 <- YTCNhatFull[j,]
          wwA <- which(regJ!=0)
          pos <- abs(zz1-(cMean))>tt
          TP <- sum(pos[wwA], na.rm=TRUE)
          ww <- which(regJ==0)     
          FP <- sum(pos[ww], na.rm=TRUE)
        }       
        return(c(FP=FP,TP=TP, pos=length(wwA), neg=length(ww)))
      }), na.rm=T)
      se <- getTPTN["TP"]/getTPTN["pos"]
      sp <- getTPTN["FP"]/getTPTN["neg"]
      return(list(tp=se,fp=sp))
    })
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
ggsave(gAUC, filename=file.path(pathFig, "AUCs.pdf"), width=10, height=7)



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

