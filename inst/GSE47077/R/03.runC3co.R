##########################################################################
# Run c3co
##########################################################################
## Run 00.preprocessing.R, 01.CallGenotypes.R, 03.TumorBoost.R before running this file
library("c3co")
patient <- "RK29"
path <- "c3co-data"
path <- R.utils::Arguments$getReadablePath(path)
filename <- sprintf("dat-%s.rds", patientID)
pathname <- file.path(path, filename)
datList <- readRDS(pathname)

dat <- readRDS(sprintf("%s/dat-%s.rds", path ,patient))
lambda.grid <- c(2e-6,1e-5,2e-5) ## penalty
lambda.grid <- seq(from=1e-6, to=1e-4, length=10)
p.list <- 2:length(datList) ## candidate number of subclones
parameters.grid <- list(lambda1=lambda.grid, lambda2=lambda.grid, nb.arch = p.list)
resC3co <- c3co(datList, pathSeg="segDat.rds", parameters.grid = parameters.grid, verbose=TRUE)
