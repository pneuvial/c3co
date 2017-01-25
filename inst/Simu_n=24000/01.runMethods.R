message("Run C3CO and FLLat on all data sets created by 00.setup.R")
mc.cores <- 1
message(sprintf("Parallelization on %s", mc.cores))

deconv <- function(b, stat) {
  message(sprintf("Run C3CO on %s, for data set %s", stat, b))
  pathArch <- Arguments$getWritablePath(sprintf("archetypeData%s_%s_c3co", stat, framework))
  pathc3co <- Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
  filename <- sprintf("archData_B=%s_%s.rds", b, "c3co")
  filepath <- file.path(pathc3co,filename)
  if(!file.exists(filepath)|| forcec3co){
    dat <- readRDS(sprintf("%s/dat_B=%s.rds",pathSim,b))
    casRes <- c3co(dat, nb.arch.grid=p.list, stat=stat,output.dir= pathc3co)
    saveRDS(casRes, filepath)
  }else{
  message("file already exists")
  }
}

c <- mclapply(B, deconv, stat="C1C2", mc.cores = mc.cores)
c <- mclapply(B, deconv, stat="TCN", mc.cores = mc.cores)

pathArch <- Arguments$getWritablePath(sprintf("archetypeData%s_%s_FLLAT", "TCN", framework))
force=FALSE
deconvFLLAT <- function(b) {
  message(sprintf("Run FLLat on %s, for data set %s", "TCN", b))
  pathSim <- Arguments$getWritablePath(sprintf("simData"))
  dat <- readRDS(sprintf("%s/dat_B=%s.rds",pathSim,b))
  YTCNtoSeg <- t(sapply(dat, function(cc) cc$tcn))
  message(sprintf("TCN is transformed to log2"))
  Y <- log2(YTCNtoSeg)-1
  pathfllat <- Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
  filename <- sprintf("archData_B=%s_%s.rds", b, "FLLAT")
  filepath <- file.path(pathfllat,filename)
  if(!file.exists(filepath)|| force){
    result.pve <- FLLat.PVE(na.omit(t(Y)),J.seq=p.list) 
    result.FLLAT <- sapply(p.list, function(pp) {
      rr <- FLLat.BIC(na.omit(t(Y)),J=pp)
      W <- rr$opt.FLLat$Theta
      Z <- rr$opt.FLLat$Beta   
      resFLLAT <- new("c3coClass", BIC=rr$opt$bic, PVE=result.pve$PVEs[which(pp==p.list)], res=new("posFused", S=list(Z=Z), W=t(W),E=list(Y=t(W)%*%t(Z))), param=list(nb.arch=pp), bkp=list(NULL))
      saveRDS(resFLLAT,sprintf("%s/featureData,p=%s.rds",pathfllat,pp))
      return(resFLLAT)
    })
    saveRDS(result.FLLAT, filepath)
  }else{
    message("file already exists")
  }
}
c <- mclapply(B, deconvFLLAT, mc.cores = mc.cores)


