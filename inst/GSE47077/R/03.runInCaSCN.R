##########################################################################
# Run InCaSCN
##########################################################################
## Run 00.preprocessing.R, 01.CallGenotypes.R, 03.TumorBoost.R before running this file
library(InCaSCN)
patient <- "RK29"
path <- "data"
dat <- readRDS(sprintf("%s/dat-%s.rds", path ,patient))
output.dir <- Arguments$getWritablePath(sprintf("resultsInCaSCN-%s",patient))

resInCaSCN <- InCaSCN(dat,lambda1.grid=c(2e-6,1e-5,2e-5),lambda2.grid=c(2e-6,1e-5,2e-5), output.dir=output.dir,nb.arch.grid = 2:9)


