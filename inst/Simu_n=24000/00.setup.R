### Package to load 
library(c3co)
library(FLLat)
library(parallel)
library(R.utils)
library(acnr)
## Number of subclones
p.list <- 2:15
## Type of simulations
framework <- "realistic"
## Force simulation and weights
forceM <- forceSim <- FALSE
## Number of simulated data sets
Bmax <- 100

### Parameters
len <- 800*3
## 3 is to obtain around 8000 point for heterozygous
nBkp <- 10

message("Path to save subclones : simArchData")
pathSimArch <- Arguments$getWritablePath("simArchData")
filename <- file.path(pathSimArch, "subClones.rds")
nbClones <- 5

####################################################################
## parameters of subclones and samples
####################################################################
set.seed(10)


if (!file.exists(filename)){
    message("subClones.rds doesn't exist in simArchData directory")
    subClones <- createSubclones(nbClones, len, nBkp)
    saveRDS(subClones, filename)
  }else {
    subClones <- readRDS(filename)
    message("subClones.rds has been loaded")
}

### Simulate Samples
message("Simulations of subclones")
n <- 30
forcec3co <- TRUE

pathWeight <- Arguments$getWritablePath(sprintf("weightData"))
B <- 1:Bmax
for (bb in B) {
  filename <- sprintf("weight,n=%s,b=%s.rds", n,bb)
  runWeights <- (forceM || !file.exists(file.path(pathWeight, filename)))
  if (runWeights) {
    message(sprintf("Create and save weight matrix for sample %s", bb))
    M <- getWeightMatrix(70,20, nbClones, n)
    saveRDS(M,file = file.path(pathWeight, filename))
  }else{
    message(sprintf("Load weight matrix for sample %s", bb))
    M <- readRDS(file.path(pathWeight, filename)) 
  }
  pathSim <- Arguments$getWritablePath(sprintf("simData"))
  fileSim <- sprintf("%s/dat_B=%s.rds",pathSim,bb)
  runSim <- (forceSim || !file.exists(fileSim))
  if (runSim) {    
### simulations of samples
    message(sprintf("Create and save DNA copy number profile for sample %s", bb))
    dat <- apply(M, 1, mixSubclones, subClones = subClones, fracN = NULL)
    saveRDS(dat, fileSim)
  }else{
    message(sprintf("Load DNA copy number profile for sample %s", bb))
    dat <- readRDS(sprintf("%s/dat_B=%s.rds",pathSim,bb))
  }
}
