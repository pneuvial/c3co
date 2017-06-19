library("c3co")
library("FLLat")
library("R.utils")
library("acnr")
library("future")
library("listenv")

plan(multiprocess, workers=3L) ## see https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html

lambda.grid <- 10^(-seq(from=6, to=4, length.out=10))
#lambda.grid <- seq(from=1e-5, to=1e-1, length=10)
p.list <- 2:10 ## candidate number of subclones
parameters.grid <- list(lambda=lambda.grid, nb.arch=p.list)
forceM <- forceSim <- FALSE

set.seed(10)  ## for reproducibility

####################################################################
## parameters of subclones and samples
####################################################################
n <- 30
nbClones <- 5
nbSimu <- 100

len <- 8000*3  ## to obtain ~800 heterozygous loci
nbClones <- 5
nBkp <- 10     ## Breakpoints in subclones

####################################################################
## paths to save results
####################################################################
pathRes <- Arguments$getWritablePath("results")
pathWeights <- Arguments$getWritablePath("weights")
pathSubClones <- Arguments$getWritablePath("subclones")
pathDat <- Arguments$getWritablePath("data")
