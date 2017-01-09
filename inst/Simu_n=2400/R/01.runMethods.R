deconv <- function(b, stat) {
### Weights
  pathArch <- Arguments$getWritablePath(sprintf("archetypeData%s_%s_c3co", stat, framework))
  pathc3co <- Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
  filename <- sprintf("archData_B=%s_%s.rds", b, "c3co")
  filepath <- file.path(pathc3co,filename)
  if(!file.exists(filepath)|| forcec3co){
    casRes <- c3co(dat, nb.arch.grid=p.list, stat=stat,output.dir= pathc3co)
    saveRDS(casRes, filepath)
  }
}
c <- mclapply(B, deconv, stat="C1C2", mc.cores = 20)
c <- mclapply(B, deconv, stat="TCN", mc.cores = 20)
library(FLLat)
pathArch <- Arguments$getWritablePath(sprintf("archetypeData%s_%s_FLLAT", "TCN", framework))
deconvFLLAT <- function(b) {
### Weights
  pathSim <- Arguments$getWritablePath(sprintf("simData"))
  dat <- readRDS(sprintf("%s/dat_B=%s.rds",pathSim,b))
  YTCNtoSeg <- t(sapply(dat, function(cc) cc$tcn))
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
      resFLLAT <- list(res=list(Z=Z, W=t(W),Y.hat=(Y=t(W)%*%t(Z))), PVE=result.pve[which(pp==p.list)])
      
      saveRDS(resFLLAT,sprintf("%s/featureData,p=%s.rds",pathfllat,pp))
    })
    saveRDS(result.FLLAT, filepath)
  }else{
    print("file already exists")
  }
}
c <- mclapply(B, deconvFLLAT, mc.cores = 20)


