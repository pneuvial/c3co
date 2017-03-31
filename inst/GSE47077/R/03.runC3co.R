##########################################################################
# Run c3co
##########################################################################
## Run 00.preprocessing.R, 01.CallGenotypes.R, 03.TumorBoost.R before running this file
library("c3co")
path <- R.utils::Arguments$getReadablePath("c3co-data")
filename <- sprintf("dat-%s.rds", patientID)
pathname <- file.path(path, filename)
datList <- readRDS(pathname)

lambda.grid <- seq(from=1e-6, to=1e-4, length=10)
p.list <- 2:length(datList) ## candidate number of subclones
parameters.grid <- list(lambda1=lambda.grid, lambda2=lambda.grid, nb.arch = p.list)

resC3co <- c3co(datList, pathSeg="segDat.rds", parameters.grid = parameters.grid, verbose=TRUE)
#seg <- resC3co@segDat
#seg$bkp <- resC3co@bkp
#saveRDS(seg, "segDat.rds")
