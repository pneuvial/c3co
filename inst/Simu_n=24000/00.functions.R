############################################################
## Some functions to load data and compute ROC and AUCs 
############################################################
createSubclones <- function(nbClones, len, nBkp) {
  ### Breakpoints in subclones
  interval <- 1:(len - 1)
  u <- numeric(0)
  minLength <- 10
  while (length(u) < nBkp) {
    j <- sample(x = interval, size = 1, replace = FALSE)
    end <- j
    u <- c(u, j)
    b.inf <- max(1, j - minLength)
    b.sup <- min(len, j + minLength)
    v <- b.inf:b.sup
    interval <- setdiff(interval, v)
  }
  
  e <- c(sort(sample(size = nbClones-1, x = 1:(length(u) - 1), replace = FALSE)), length(u))
  s <- c(1,e + 1)[-(length(e) + 1)]
  message("Breakpoints by clones")
  bkpsByClones <- mapply(function(ss,ee){
    sort(u[ss:ee])
  }, s,e)
  o <- order(u)
  orderBkpsByClones <- mapply(function(ss,ee){
    sort(o[ss:ee])
  }, s,e)
  ### Regions
  
  regNames <- c("(1,1)","(0,1)","(1,2)","(0,2)")
  pattern <- "\\(([0-9]),([0-9])\\)"
  regAnnot <- data.frame(region = regNames, 
                         freq = rep(1/4,4), stringsAsFactors = FALSE)
  regAnnot$C1 <- as.numeric(gsub(pattern, "\\1", regNames))
  regAnnot$C2 <- as.numeric(gsub(pattern, "\\2", regNames))
  message("Regions by subclones")
  candidateRegions <- function(regName) {
    if (is.null(regName)) 
      return(regAnnot[, "region"])
    reg <- subset(regAnnot, region == regName)
    d1 <- regAnnot[, "C1"] - reg[, "C1"]
    d2 <- regAnnot[, "C2"] - reg[, "C2"]
    ww <- which((d1 & !d2) | (!d1 & d2))
    regAnnot[ww, "region"]
  }
  
  tmp <- matrix("(1,1)",nrow = nbClones, ncol = len)
  
  for (rr in 1:nrow(tmp)) {
    start <- c(1,bkpsByClones[[rr]] + 1)
    end <- c(bkpsByClones[[rr]], len)
    for (bb in 1:length(start)) {
      if(bb == 1) {
        reg = NULL
      }else {
        reg <- tmp[rr,start[bb - 1]]
      }      
      candReg <- candidateRegions(reg)
      reg <- sample(size = 1, candReg)
      tmp[rr,start[bb]:end[bb]] <- reg
    }
  }
  
  regionsByClones <- sapply(1:length(bkpsByClones), function(bb){
    bkp <- c(bkpsByClones[[bb]],len)
    tmp[bb,bkp]
  })
  
  c1c2 <- "\\(([0-9]),([0-9])\\)"
  dataAnnotTP <- acnr::loadCnRegionData(dataSet = "GSE13372", tumorFrac = 1)
  dataAnnotN <- acnr::loadCnRegionData(dataSet = "GSE13372", tumorFrac = 0)
  message("Build subclones")
  subClones <- buildSubclones(len, dataAnnotTP, dataAnnotN,
                              nbClones, bkpsByClones,regionsByClones)
  return(subClones)
  
}

runFLLAT <- function(dat, p.list) {
  message(sprintf("TCN is transformed to log2"))
  YTCNtoSeg <- t(sapply(dat, function(cc) cc$tcn))
  Y <- log2(YTCNtoSeg) - 1
  result.pve <- FLLat::FLLat.PVE(na.omit(t(Y)), J.seq = p.list) 
  result.FLLAT <- sapply(p.list, function(pp) {
    rr <- FLLat::FLLat.BIC(na.omit(t(Y)), J = pp)
    W <- rr$opt.FLLat$Theta
    Z <- rr$opt.FLLat$Beta   
    resFLLAT <- new("c3coClass", BIC = rr$opt$bic, PVE = result.pve$PVEs[which(pp == p.list)], res = new("posFused", S = list(Z = Z), W = t(W),E = list(Y = t(W) %*% t(Z))), param = list(nb.arch = pp), bkp = list(NULL))
    #saveRDS(resFLLAT,sprintf("%s/featureData,p=%s.rds",pathfllat,pp))
    return(resFLLAT)
  })
  result.FLLAT
}


library(R.utils)
loadDataBest <- function(mm, stat, framework, b){
  print(sprintf("meth=%s, stat=%s", mm, stat))
  pathArch <- Arguments$getWritablePath(sprintf("archetypeData%s_%s_%s", stat, framework, mm))
  pathMeth <- Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
  filename <- sprintf("archData_B=%s_%s.rds", b, mm)
  print(filename)
  file <- file.path(pathMeth,filename)
  dataBest <- NULL
  if (file.exists(file)) {
    res <- readRDS(file)
      pves <- unlist(sapply(res, function(rr) rr@PVE))
      pBest <- min(c(which(diff(pves) < 1e-3),length(pves) ))
      dataBest <- res[[pBest]]
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


SESPC1C2 <- function(Z1,Z2,alteredLoci,ind, tol, c1Mean, c2Mean){
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
      return(c(FP = FP,TP = TP, pos = length(wwA), neg = length(ww)))
    }))
    se <- getTPTN["TP"]/getTPTN["pos"]
    sp <- getTPTN["FP"]/getTPTN["neg"]
    return(list(tp = se,fp = sp))
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
      
      return(c(FP = FP,TP = TP, pos = length(wwA), neg = length(ww)))
    }), na.rm = T)
    se <- getTPTN["TP"]/getTPTN["pos"]
    sp <- getTPTN["FP"]/getTPTN["neg"]
    return(list(tp = se,fp = sp))
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



getPos <- function(tt, zz1,cMean){
  pos <- (abs(zz1 - cMean) >= tt)
  return(pos)
}


getPosC1C2 <- function(tt,zz1,zz2, c1Mean, c2Mean){
  pos <- ((abs(zz1 - c1Mean) >= tt | abs(zz2 - c2Mean) >= tt))
  return(pos)
}
