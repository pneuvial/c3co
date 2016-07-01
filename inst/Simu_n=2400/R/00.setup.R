library(InCaSCN)
library(FLLat)
library(parallel)
## Number of archetypes
p.list <- 2:15
framework <- "realistic"
stat <- "C1C2"
forceM <- forceSim <- FALSE


####################################################################
  ## parameters of subclones and samples
####################################################################

pathSimArch <- Arguments$getWritablePath("simArchData")
filename <- file.path(pathSimArch, "subClones.rds")

if(!file.exists(filename)){
  set.seed(10)
### Parameters
  len <- 800*3
  ## 3 is to obtain around 8000 point for heterozygous
  nbClones <- 5
  nBkp <- 10
### Breakpoints in subclones
  interval <- 1:(len - 1)
  u <- numeric(0)
  minLength <- 100
  while (length(u) < nBkp) {
    j <- sample(x = interval, size = 1, replace = FALSE)
    jend <- j
    u <- c(u, j)
    b.inf <- max(1, j - minLength)
    b.sup <- min(len, j + minLength)
    v <- b.inf:b.sup
    interval <- setdiff(interval, v)
  }

  e <- c(sort(sample(size=4, x=1:(length(u)-1), replace=FALSE)), length(u))
  s <- c(1,e+1)[-(length(e)+1)]
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
  candidateRegions <- function(regName) {
    if (is.null(regName)) 
      return(regAnnot[, "region"])
    reg <- subset(regAnnot, region == regName)
    d1 <- regAnnot[, "C1"] - reg[, "C1"]
    d2 <- regAnnot[, "C2"] - reg[, "C2"]
    ww <- which((d1 & !d2) | (!d1 & d2))
    regAnnot[ww, "region"]
  }

  tmp <- matrix("(1,1)",nrow=nbClones, ncol=len)

  for(rr in 1:nrow(tmp)){
    start <- c(1,bkpsByClones[[rr]]+1)
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

  c1c2 <- "\\(([0-9]),([0-9])\\)"
  dataAnnotTP <- loadCnRegionData(dataSet="GSE13372", tumorFrac=1)
  dataAnnotN <- loadCnRegionData(dataSet="GSE13372", tumorFrac=0)
  subClones <- buildSubclones(len,dataAnnotTP, dataAnnotN,
                              nbClones, bkpsByClones,regionsByClones)
  saveRDS(subClones, filename)
}else{
  subClones <- readRDS(filename)
}

### Simulate Samples
n <- 30
forceInCaSCN <- FALSE

pathWeight <- Arguments$getWritablePath(sprintf("weightData"))
B <- 1:100
fir(bb in B){
  filename <- sprintf("weight,n=%s,b=%s.rds", n,b)
  runWeights <- (forceM||!file.exists(file.path(pathWeight, filename)))
  if(runWeights){
    M <- getWeightMatrix(70,20, nbClones, n)
    saveRDS(M,file=file.path(pathWeight, filename))
  }else{
    M <- readRDS(file.path(pathWeight, filename)) 
  }
  pathSim <- Arguments$getWritablePath(sprintf("simData"))
  fileSim <- sprintf("%s/dat_B=%s.rds",pathSim,b)
  runSim <- (forceSim||!file.exists(fileSim))
  if(runSim){    
### simulations of samples
    dat <- apply(M, 1, mixSubclones, subClones=subClones, fracN=NULL)
    saveRDS(dat, fileSim)
  }else{
    dat <- readRDS(sprintf("%s/dat_B=%s.rds",pathSim,b))
  }
}
