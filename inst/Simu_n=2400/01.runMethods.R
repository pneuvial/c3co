deconv <- function(b, stat, nb.arch.grid) {
### Weights
  pathArch <- Arguments$getWritablePath(sprintf("archetypeData%s_%s_c3co", stat, framework))
  pathc3co <- Arguments$getWritablePath(sprintf("%s/features_B=%s/", pathArch, b))
  filename <- sprintf("archData_B=%s_%s.rds", b, "c3co")
  filepath <- file.path(pathc3co,filename)
  if (!file.exists(filepath)|| forcec3co){
    dat <- readRDS(sprintf("%s/dat_B=%s.rds", pathSim, bb))
    casRes <- c3co(dat, nb.arch.grid=nb.arch.grid, stat=stat, output.dir=pathc3co)
    saveRDS(casRes, filepath)
  }
}
dummy <- mclapply(B, deconv, stat="C1C2", nb.arch.grid=p.list, mc.cores = mc.cores)
dummy <- mclapply(B, deconv, stat="TCN", nb.arch.grid=p.list, mc.cores = mc.cores)

pathArch <- Arguments$getWritablePath(sprintf("archetypeData%s_%s_FLLAT", "TCN", framework))
deconvFLLAT <- function(b) {
### Weights
  pathSim <- Arguments$getWritablePath("simData")
  dat <- readRDS(sprintf("%s/dat_B=%s.rds", pathSim, b))
  YTCNtoSeg <- t(sapply(dat, function(cc) cc$tcn))
  Y <- log2(YTCNtoSeg)-1
  pathfllat <- Arguments$getWritablePath(sprintf("%s/features_B=%s/", pathArch, b))
  filename <- sprintf("archData_B=%s_%s.rds", b, "FLLAT")
  filepath <- file.path(pathfllat,filename)
  if(!file.exists(filepath)|| force){
    result.pve <- FLLat.PVE(na.omit(t(Y)),J.seq=p.list) 
    result.FLLAT <- sapply(p.list, function(pp) {
      rr <- FLLat.BIC(na.omit(t(Y)),J=pp)
      W <- rr$opt.FLLat$Theta
      Z <- rr$opt.FLLat$Beta   
      resFLLAT <- list(res=list(Z=Z, W=t(W),Y.hat=(Y=t(W)%*%t(Z))), PVE=result.pve[which(pp==p.list)])
      
      saveRDS(resFLLAT,sprintf("%s/featureData,p=%s.rds",pathfllat,pp))
    })
    saveRDS(result.FLLAT, filepath)
  }else{
    print("file already exists")
  }
}
dummy <- mclapply(B, deconvFLLAT, mc.cores = mc.cores)


