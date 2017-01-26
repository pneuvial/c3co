library("c3co")
library("FLLat")
library("parallel")
library("R.utils")
library("acnr")


p.list <- 2:15  ## candidate number of subclones
forceM <- forceSim <- FALSE
mc.cores <- 20L

set.seed(10)  ## for reproducibility

####################################################################
## parameters of subclones and samples
####################################################################
n <- 30
nbClones <- 5
nbSimu <- 10

len <- 800*3  ## 3 is to obtain around 800 point for heterozygous
nbClones <- 5
nBkp <- 10     ## Breakpoints in subclones

