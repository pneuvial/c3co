library(InCaSCN)
patient <- "OV03-17"

output.dir <- Arguments$getWritablePath(sprintf("resultsInCaSCN-%s",patient))
resInCaSCN <- InCaSCN(NULL,lambda1.grid=c(1e-6,2e-6,1e-5,2e-5),lambda2.grid=c(1e-6,2e-6,1e-5,2e-5), output.dir=output.dir,nb.arch.grid = 2:7)

