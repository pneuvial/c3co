### Mesure Performance 
library("c3co")
library("FLLat")
library("R.utils")
library("acnr")
library("future")
library("listenv")

plan(multiprocess, workers=3L) ## see https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html

lambda.grid <- seq(from=1e-6, to=1e-4, length=10)
p.list <- 2:10 ## candidate number of subclones
parameters.grid1 <- list(lambda=lambda.grid, nb.arch=p.list)
parameters.grid2 <- list(lambda1=lambda.grid,lambda2=lambda.grid, nb.arch=p.list)
forceM <- forceSim <- FALSE

set.seed(10)  ## for reproducibility

####################################################################
## parameters of subclones and samples
####################################################################
n <- 30
nbClones <- 5
nbSimu <- 10

len <- 800*3  ## to obtain ~800 heterozygous loci
nbClones <- 5
nBkp <- 10     ## Breakpoints in subclones

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


subClones <- simulateSubclones(len, nbClones, nBkp)

## simulate copy number profiles from subclones and weight matrices
weightMats <- list()
dats <- list()
resC1C2_1 <- listenv::listenv()
resC1C2_2 <- listenv::listenv()
for (ss in 1:nbSimu) {
  ## weight matrix
  M <- rSparseWeightMatrix(n, nbClones, 0.7)
  weightMats[[ss]] <- M
  
  ## simulated profiles
  dat <- mixSubclones(subClones=subClones, M)
  dats[[ss]] <- dat
  ## c3co
  resC1C2_1[[ss]] %<-% c3co(dat, parameters.grid=parameters.grid1, stat="C1C2")
  resC1C2_2[[ss]] %<-% c3co(dat, parameters.grid=parameters.grid2, stat="C1C2")
}

### Comparison
list_res <- as.list(resC1C2_2)
nbEq <- sum(sapply(list_res, function (rr){
  res <- sapply(rr@fit, function(ff){
    r <- ff@param$lambda1==ff@param$lambda2
  })
}))/(length(p.list)*nbSimu)

## Finally the best results with BIC are obtained when lambda1=lambda2
## To save time we can use in experiments lambda1=lambda2.



