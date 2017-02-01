library("c3co")
library("FLLat")
library("parallel")
library("R.utils")
library("acnr")
library("future")
library("listenv")

plan(multiprocess, workers=3L) ## see https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html

p.list <- 2:10 ## candidate number of subclones
forceM <- forceSim <- FALSE

set.seed(10)  ## for reproducibility

####################################################################
## parameters of subclones and samples
####################################################################
n <- 30
nbClones <- 5
nbSimu <- 5

len <- 800*3  ## 3 is to obtain around 800 point for heterozygous
nbClones <- 5
nBkp <- 10     ## Breakpoints in subclones

####################################################################
## paths to save results
####################################################################
pathRes <- Arguments$getWritablePath("results")
pathWeights <- Arguments$getWritablePath("weights")
pathSubClones <- Arguments$getWritablePath("subclones")
pathDat <- Arguments$getWritablePath("data")


