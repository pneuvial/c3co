##########################################################################
# Run c3co
##########################################################################
## Run 00.preprocessing.R, 01.CallGenotypes.R, 03.TumorBoost.R before running this file
library("c3co")
patient <- "RK29"
path <- "data"
dat <- readRDS(sprintf("%s/dat-%s.rds", path ,patient))
output.dir <- sprintf("results_c3co-%s",patient)

grid <- c(2e-6,1e-5,2e-5)
resC3co <- c3co(dat,lambda1.grid=grid, lambda2.grid=grid, 
                saveResults=TRUE, output.dir=output.dir,nb.arch.grid = 2:9)


