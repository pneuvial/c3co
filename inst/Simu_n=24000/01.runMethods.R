message("Run C3CO and FLLat on all data sets created by 00.setup.R")
mc.cores <- 1L
message(sprintf("Parallelization on %s", mc.cores))

deconv <- function(b) {
  message(sprintf("Run C3CO on %s, for data set %s", "C1C2", b))
  pathArch <- R.utils::Arguments$getWritablePath(sprintf("archetypeData%s_%s_c3co", "C1C2", framework))
  pathc3co <- R.utils::Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
  filename <- sprintf("archData_B=%s_%s.rds", b, "c3co")
  filepath <- file.path(pathc3co,filename)
  if (!file.exists(filepath) || forcec3co) {
    casRes <- c3co(readRDS(sprintf("%s/dat_B=%s.rds",pathSim,b)), nb.arch.grid = p.list, stat = "C1C2" ,output.dir = pathc3co)
    saveRDS(casRes, filepath)
  }else{
    message("file already exists")
  }
  message(sprintf("Run C3CO on %s, for data set %s", "TCN", b))
  pathArch <- R.utils::Arguments$getWritablePath(sprintf("archetypeData%s_%s_c3co", "TCN", framework))
  pathc3co <- R.utils::Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
  filename <- sprintf("archData_B=%s_%s.rds", b, "c3co")
  filepath <- file.path(pathc3co,filename)
  if (!file.exists(filepath) || forcec3co) {
    casRes <- c3co(readRDS(sprintf("%s/dat_B=%s.rds",pathSim,b)), nb.arch.grid = p.list, stat = "TCN" ,output.dir = pathc3co)
    saveRDS(casRes, filepath)
    
  }else{
    message("file already exists")
  }
  
}


pathArch <- R.utils::Arguments$getWritablePath(sprintf("archetypeData%s_%s_FLLAT", "TCN", framework))
force = FALSE
deconvFLLAT <- function(b) {
  message(sprintf("Run FLLat on %s, for data set %s", "TCN", b))
  pathSim <- R.utils::Arguments$getWritablePath(sprintf("simData"))
  pathfllat <- R.utils::Arguments$getWritablePath(sprintf("%s/features_B=%s/",pathArch,b))
  filename <- sprintf("archData_B=%s_%s.rds", b, "FLLAT")
  pathFeat <- file.path(pathfllat,filename)
  dat <- readRDS(sprintf("%s/dat_B=%s.rds",pathSim,b))
  if (!file.exists(pathFeat) || force) {
    res <- runFLLAT(dat, p.list)
    saveRDS(res, pathFeat)
  }else{
    message("file already exists")
  }
}


dum <- parallel::mclapply(B, deconv, mc.cores = mc.cores)
dum <- parallel::mclapply(B, deconvFLLAT, mc.cores = mc.cores)


