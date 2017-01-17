library(R.utils)
loadDataBest <- function(mm, stat, framework, b){
  print(sprintf("meth=%s, stat=%s", mm, stat))
  pathArch <- Arguments$getWritablePath(sprintf("archetypeData%s_%s_%s", stat, framework, mm))
  pathMeth <- Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
  filename <- sprintf("archData_B=%s_%s.rds", b, mm)
  print(filename)
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
      pBest <- min(c(which(diff(listRes$PVE)<1e-3),length(listRes$PVE) ))
      dataBest <- listRes$res[[pBest]]
    }
    if(mm=="c3co"){
      listRes <- list()
      listRes$PVE <- unlist(sapply(res, function (rr) rr$PVE))
      listRes$nb.arch <- unlist(sapply(res, function (rr) rr$param$nb.arch))
      listRes$res <- lapply(res, function (rr) rr$res)
      pBest <- min(c(which(diff(listRes$PVE)<1e-3),length(listRes$PVE) ))
      dataBest <- listRes$res[[pBest]]
      dataBest$bkp <- res[[pBest]]$bkp[[1]]
    }
    
  }
  return(dataBest)
}
expand <- function(mat, start, end){
  matEx <- t(apply(mat, 2, function(cc){
    y <- unlist(sapply(seq(along=start), function(ii){
      rep(cc[ii], times=(end[ii]-start[ii]+1))
    }))
  }))
  return(matEx)
}


ComputeROC <- function(roc, FPRs) {
  TPRs <- sapply(FPRs, FUN=function(fp) {
      ww <- which(roc["fp",]<=fp)
      max(roc["tp",ww])
  })
  TPRs[is.infinite(TPRs)] <- 0
   return(TPRs)
}


SESPC1C2 <- function(Z1,Z2,alteredLoci,ind, tol, c1Mean, c2Mean){
  sespRes <- sapply(tol, function(tt){
    getTPTN <- rowSums(sapply(1:ncol(alteredLoci), function(j){
      k <- ind[j]
      regJ <- alteredLoci[,j]
      zz1 <- Z1[k,]
      zz2 <- Z2[k,]
      wwA <- which(regJ)
      pos <- ((abs(zz1-c1Mean)>=tt | abs(zz2-c2Mean)>=tt))
      TP <- sum(pos[wwA])
      ww <- which(!regJ)      
      FP <- sum(pos[ww])  
      return(c(FP=FP,TP=TP, pos=length(wwA), neg=length(ww)))
    }))
    se <- getTPTN["TP"]/getTPTN["pos"]
    sp <- getTPTN["FP"]/getTPTN["neg"]
    return(list(tp=se,fp=sp))
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
      pos <- (abs(zz1-cMean)>=tt)
      TP <- sum(pos[wwA])
      ww <- which(!regJ)      
      FP <- sum(pos[ww])
      
      return(c(FP=FP,TP=TP, pos=length(wwA), neg=length(ww)))
    }), na.rm=T)
    se <- getTPTN["TP"]/getTPTN["pos"]
    sp <- getTPTN["FP"]/getTPTN["neg"]
    return(list(tp=se,fp=sp))
  })
  return(sespRes)
}


ComputeAUC <- function(roc, FPRs) {
    
  TPRs <- sapply(FPRs, FUN=function(fp) {
    ww <- which(roc["fp",]<=fp)
    max(roc["tp",ww])
  })
  TPRs[is.infinite(TPRs)] <- NaN
  auc <- sum(lintegrate(FPRs, TPRs, xint=FPRs)) 
  return(auc)
}



getPos <- function(tt, zz1,cMean){
  pos <- (abs(zz1-cMean)>=tt)
  return(pos)
}


getPosC1C2 <- function(tt,zz1,zz2, c1Mean, c2Mean){
  pos <- ((abs(zz1-c1Mean)>=tt | abs(zz2-c2Mean)>=tt))
  return(pos)
}
