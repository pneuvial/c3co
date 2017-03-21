##########################################################################
# Run c3co
##########################################################################
## Run 00.preprocessing.R, 01.CallGenotypes.R, 03.TumorBoost.R before running this file
library("c3co")
patient <- "RK29"
path <- "data"
dat <- readRDS(sprintf("%s/dat-%s.rds", path ,patient))

lambda.grid <- c(2e-6,1e-5,2e-5) ## penalty
p.list <- 2:9 ## candidate number of subclones
parameters.grid <- list(lambda1=lambda.grid, lambda2=lambda.grid, nb.arch = p.list)
resC3co <- c3co(dat, parameters.grid = parameters.grid, verbose=TRUE)
