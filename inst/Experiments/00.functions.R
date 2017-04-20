simulateSubclones <- function(len, nbClones, nBkp) {
    interval <- 1:(len - 1)
    u <- numeric(0)
    minLength <- 100
    while (length(u) < nBkp) {
        j <- sample(x=interval, size=1, replace=FALSE)
        u <- c(u, j)
        b.inf <- max(1, j - minLength)
        b.sup <- min(len, j + minLength)
        v <- b.inf:b.sup
        interval <- setdiff(interval, v)
    }
    
    e <- c(sort(sample(size=nbClones-1, x=1:(length(u)-1), replace=FALSE)),length(u))
    s <- c(1,e+1)[-(length(e)+1)]
    bkpsByClones <- mapply(function(ss,ee){
        sort(u[ss:ee])
    }, s, e)
    o <- order(u)
    
    ### Regions
    regNames <- c("(1,1)","(0,1)","(1,2)","(0,2)")
    pattern <- "\\(([0-9]),([0-9])\\)"
    regAnnot <- data.frame(region = regNames, 
                           freq = rep(1/4,4), stringsAsFactors = FALSE)
    regAnnot$C1 <- as.numeric(gsub(pattern, "\\1", regNames))
    regAnnot$C2 <- as.numeric(gsub(pattern, "\\2", regNames))
    
    candidateRegions <- function(regName) {
        if (is.null(regName)) 
            return(regAnnot[, "region"])
        reg <- subset(regAnnot, region == regName)
        d1 <- regAnnot[, "C1"] - reg[, "C1"]
        d2 <- regAnnot[, "C2"] - reg[, "C2"]
        ww <- which((d1 & !d2) | (!d1 & d2))
        regAnnot[ww, "region"]
    }
    
    tmp <- matrix("(1,1)", nrow=nbClones, ncol=len)
    
    for (rr in 1:nrow(tmp)) {
        start <- c(1, bkpsByClones[[rr]]+1)
        end <- c(bkpsByClones[[rr]], len)
        for(bb in 1:length(start)){
            if(bb==1){
                reg=NULL
            }else{
                reg <- tmp[rr,start[bb-1]]
            }      
            candReg <- candidateRegions(reg)
            reg <- sample(size=1, candReg)
            tmp[rr,start[bb]:end[bb]] <- reg
        }
    }
    
    regionsByClones <- sapply(1:length(bkpsByClones), function(bb){
        bkp <- c(bkpsByClones[[bb]],len)
        tmp[bb,bkp]
    })
    
    dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE13372", tumorFraction=1)
    dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE13372", tumorFraction=0)
    subClones <- c3co::buildSubclones(len, dataAnnotTP, dataAnnotN,
                                nbClones, bkpsByClones,regionsByClones)
    subClones
}


fllat <- function(dat, nb.arch.grid) {
    YTCNtoSeg <- t(sapply(dat, function(cc) cc$tcn))
    Y <- log2(YTCNtoSeg)-1
    result.pve <- FLLat::FLLat.PVE(na.omit(t(Y)), J.seq=nb.arch.grid) 
    res <- sapply(nb.arch.grid, function(pp) {
        rr <- FLLat::FLLat.BIC(na.omit(t(Y)), J=pp)
        W <- rr$opt.FLLat$Theta
        Z <- rr$opt.FLLat$Beta   
        resFLLAT <- methods::new("posFused", S=list(Z=Z), W=t(W), E=list(Y=t(W)%*%t(Z)), PVE=result.pve$PVEs[which(pp==nb.arch.grid)], BIC=rr$opt.FLLat$bic, param=list(nb.arch=pp, lambda1=rr$lam1, lambda2=rr$lam2))
    })
    methods::new("c3coFit", bkp=list(), segDat=list(), fit=res)
}

loadDataBest <- function(mm, stat, b, nbClones=5){
  print(sprintf("meth=%s, stat=%s", mm, stat))
  pathRes <- "results"
  filename <- sprintf("results_%s_%s.rds", mm, stat)
  print(filename)
  file <- file.path(pathRes,filename)
  dataBest <- NULL
  if (file.exists(file)) {
    res <- readRDS(file)[[b]]
    fit <- res@fit
    pves <- unlist(sapply(fit, function(rr) rr@PVE))
    pBest <- min(c(which(diff(pves) < 1e-3),length(pves) ))
    dataBest <- fit[[pBest]]
    return(dataBest)
  }
}

expand <- function(mat, start, end){
  matEx <- t(apply(mat, 2, function(cc){
    y <- unlist(sapply(seq(along = start), function(ii){
      rep(cc[ii], times = (end[ii] - start[ii] + 1))
    }))
  }))
  return(matEx)
}


ComputeROC <- function(roc, FPRs) {
  TPRs <- sapply(FPRs, FUN = function(fp) {
    ww <- which(roc["fp",] <= fp)
    max(roc["tp",ww])
  })
  TPRs[is.infinite(TPRs)] <- 0
  return(TPRs)
}


SESPC1C2 <- function(Z1, Z2, alteredLoci, ind, tol, c1Mean, c2Mean){
  sespRes <- sapply(tol, function(tt){
    getTPTN <- rowSums(sapply(1:ncol(alteredLoci), function(j){
      k <- ind[j]
      regJ <- alteredLoci[,j]
      zz1 <- Z1[k,]
      zz2 <- Z2[k,]
      wwA <- which(regJ)
      pos <- ((abs(zz1 - c1Mean) >= tt | abs(zz2 - c2Mean) >= tt))
      TP <- sum(pos[wwA])
      ww <- which(!regJ)      
      FP <- sum(pos[ww])  
      return(c(FP=FP, TP=TP, pos=length(wwA), neg=length(ww)))
    }))
    se <- getTPTN["TP"]/getTPTN["pos"]
    sp <- getTPTN["FP"]/getTPTN["neg"]
    return(list(tp=se, fp=sp))
  })
  return(sespRes)
}

SESPTCN <- function(Z,alteredLoci, ind,tol, cMean){
  sespRes <- sapply(tol, function(tt){
    getTPTN <- rowSums(sapply(1:ncol(alteredLoci), function(j){
      k <- ind[j]
      regJ <- alteredLoci[,j]
      zz1 <- Z[k,]
      wwA <- which(regJ)
      pos <- (abs(zz1 - cMean) >= tt)
      TP <- sum(pos[wwA])
      ww <- which(!regJ)      
      FP <- sum(pos[ww])
      
      return(c(FP=FP, TP=TP, pos=length(wwA), neg=length(ww)))
    }), na.rm=TRUE)
    se <- getTPTN["TP"]/getTPTN["pos"]
    sp <- getTPTN["FP"]/getTPTN["neg"]
    return(list(tp=se, fp=sp))
  })
  return(sespRes)
}


ComputeAUC <- function(roc, FPRs) {
  
  TPRs <- sapply(FPRs, FUN = function(fp) {
    ww <- which(roc["fp",] <= fp)
    max(roc["tp",ww])
  })
  TPRs[is.infinite(TPRs)] <- NaN
  auc <- sum(tis::lintegrate(FPRs, TPRs, xint = FPRs)) 
  return(auc)
}


## Get positive loci 
getPos <- function(tt, zz1, cMean){
  pos <- (abs(zz1 - cMean) >= tt)
  return(pos)
}

## Get positive loci for C1 C2
getPosC1C2 <- function(tt, zz1, zz2, c1Mean, c2Mean){
  pos <- ((abs(zz1 - c1Mean) >= tt | abs(zz2 - c2Mean) >= tt))
  return(pos)
}


## Compute Loss What-WT
lossW <- function(nbSimu, meth, stats,weightsMat){
  
  weightsArray <- array(dim = c(nbSimu, length(stats)),dimnames = list(b = 1:nbSimu, method = sprintf("%s-%s",meth,stats)))
  
  for (b in 1:nbSimu) {
    print(b)
    M <- weightsMat[[b]]
    WT <- cbind(M, 1 - rowSums(M))
    for (ii in 1:length(stats)) {
      stat = stats[ii]
      mm <- meth[ii]
      dataBest <- loadDataBest(mm, stat, b)
      if (!is.null(dataBest)) {
        W <- dataBest@W    
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
      weightsArray[b,ii] <- resLossWeights
    }
  }
  return(weightsArray)
}


## Compute randIndex(clust(What), clust(WT)
randIndW <- function(nbSimu, meth, stats, weightsMat){
  randIndexArray <- array(dim = c(nbSimu, length(stats)),dimnames = list(b = 1:nbSimu, method = sprintf("%s-%s",meth,stats)))
  for (b in 1:nbSimu) {
    M <- weightsMat[[b]]
    WT <- cbind(M, 1-rowSums(M))
    d_clust <- mclust::Mclust(t(WT), G=1:20)
    p.best <- dim(d_clust$z)[2]
    print(p.best)
    clustWT <-  cutree(hclust(dist(WT),method = "ward.D"), p.best)
    
    for (ii in 1:length(stats)) {
      stat = stats[ii]
      mm <- meth[ii]
      randIndex <- NA
      dataBest <- loadDataBest(mm, stat, b)
      if (!is.null(dataBest)) {
        clustWhat <- cutree(hclust(dist(dataBest@W), method = "ward.D"), p.best)
        randIndex <- mclust::adjustedRandIndex(clustWT, clustWhat)
      }
      randIndexArray[b,ii] <- randIndex
    }
  }
  return(randIndexArray)
}


## Get PVE values
pveEval <- function(nbSimu, meth, stats){
  PVEArray <- array(dim=c(nbSimu, length(stats), max(p.list)-1),dimnames=list(b=1:nbSimu, method=sprintf("%s-%s",meth,stats)))
  for(b in 1:nbSimu){
    print(b)
    for(ii in 1:length(stats)){
      stat=stats[ii]
      mm=meth[ii]
      pathRes <- "results"
      filename <- sprintf("results_%s_%s.rds", mm, stat)
      print(filename)
      file <- file.path(pathRes,filename)
      if(file.exists(file)){
        res <- readRDS(file)[[b]]
        PVEArray[b,ii,1:length(res@fit)] <- sapply(res@fit, function (rr) rr@PVE)
      }
    }
  }
  return (PVEArray)
}

## Compute AUC on Subclones
computeAUC <- function(nbSimu, meth, stats, tol, subClones, weightsMat, regionsByClones){
  rocArrayArchFull <- array(dim = c(nbSimu, length(stats), 2,length(tol)),dimnames = list(b = 1:nbSimu, method = sprintf("%s-%s",meth,stats), ROC = c("tp", "fp")))
  alteredLociInClones <- sapply(1:(length(regionsByClones) - 1), function(c){
    if (c != length(regionsByClones)){
      res <- subClones[[c]]$region != "(1,1)"
    }else{
      res <- rep(FALSE, len) 
    }
    return(res)
  })
  
  
  for (b in 1:nbSimu) {
    print(b)
    M <- weightsMat[[b]]
    WT <- cbind(M, 1 - rowSums(M))
    for (ss in 1:length(stats)) {
      stat <- stats[ss]
      mm <- meth[ss]      
      dataBest <- loadDataBest(mm,stat, b)
      if (!is.null(dataBest)) {
        message(sprintf("Compute ROC and AUC for method %s, var %s and data set %s" ,mm, stat, b))
        Z <- dataBest@S$Z
        W <- dataBest@W
        eps <- 0.1
        corrBestW <- apply(WT,2,function(zz){
          apply(W,2,function(zzh){
            sum(abs(zz - zzh) < eps,na.rm = TRUE)/length(zzh)
          }) 
        })
        ind <- apply(corrBestW, 2, which.max)
        
        if (stat == "C1C2") {
          Z1 <- dataBest@S$Z1
          Z2 <- dataBest@S$Z2
          pathRes <- "results"
          filename <- sprintf("results_%s_%s.rds", mm, stat)
          print(filename)
          file <- file.path(pathRes,filename)
          bkp <- readRDS(file)[[b]]@bkp[[1]]
          start <- c(ceiling(bkp)[-length(bkp)])
          end <- c(floor(bkp)[-1])
          Z1hatFull <- expand(Z1,start, end)
          Z2hatFull <- expand(Z2,start, end)
          ZhatFull <- expand(Z,start, end)
          
          SESP <- SESPC1C2(Z1hatFull,Z2hatFull,alteredLociInClones,ind, tol, 1, 1)
        }else{
          if(mm == "FLLAT") {
            ZhatFull <- t(2*2^Z)          
          }else{
            pathRes <- "results"
            filename <- sprintf("results_%s_%s.rds", mm, stat)
            print(filename)
            file <- file.path(pathRes,filename)
            bkp <- readRDS(file)[[b]]@bkp[[1]]
            start <- c(ceiling(bkp)[-length(bkp)])
            end <- c(floor(bkp)[-1])        
            ZhatFull <- expand(Z,start, end)
          }
          SESP <- SESPTCN(ZhatFull,alteredLociInClones, ind,tol, 2)
          
        }
        rocArrayArchFull[b,ss,"tp",] <- unlist(SESP["tp",])
        rocArrayArchFull[b,ss,"fp",] <- unlist(SESP["fp",])
      }
    }
  }
  FPRs <-  sort(unique(as.vector(rocArrayArchFull[,,"fp",])))
  AUCs_arch<- apply(rocArrayArchFull, 2, function(roc){
    apply(roc,1,ComputeAUC, FPRs)
  })
  
  return(AUCs_arch)
}

