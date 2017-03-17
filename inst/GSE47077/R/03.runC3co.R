##########################################################################
# Run c3co
##########################################################################
## Run 00.preprocessing.R, 01.CallGenotypes.R, 03.TumorBoost.R before running this file
library("c3co")
patient <- "RK29"
path <- "data"
dat <- readRDS(sprintf("%s/dat-%s.rds", path ,patient))
output.dir <- R.utils::Arguments$getWritablePath(sprintf("results_c3co-%s",patient))

resC3co <- c3co(dat,lambda1.grid=c(2e-6,1e-5,2e-5),lambda2.grid=c(2e-6,1e-5,2e-5), saveResults=TRUE, output.dir=output.dir,nb.arch.grid = 2:9)


