##########################################################################
# Run c3co
##########################################################################

library("c3co")
patient <- "RK29"

path <- "c3co-data"
path <- R.utils::Arguments$getReadablePath(path)
filename <- sprintf("dat-%s.rds", patientID)
pathname <- file.path(path, filename)
datList <- readRDS(pathname)

lambda.grid <- c(2e-6, 1e-5, 2e-5) ## penalty
#lambda.grid <- seq(from=1e-6, to=1e-4, length=10)
p.list <- 2:length(datList) ## candidate number of subclones
parameters.grid <- list(lambda1=lambda.grid, lambda2=lambda.grid, nb.arch = p.list)
resC3co <- c3co(datList, parameters.grid=parameters.grid, verbose=TRUE)

##resC3co <- c3co(datList, pathSeg="segDat.rds", parameters.grid = parameters.grid, verbose=TRUE)
